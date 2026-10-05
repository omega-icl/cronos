// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

// Python binding of FFDom (ocbase.hpp): the domain of an independent variable
// -- bounds, finite elements, collocation node family and count -- as consumed
// by FFModel.add_domain and add_input.
#include <pybind11/native_enum.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <sstream>

#include "ocbase.hpp"

namespace py = pybind11;

void
mc_ffdom(py::module_& m)
{
  py::class_<mc::FFDom> pyFFDom(m, "FFDom", R"doc(
Domain of an independent variable: bounds, finite elements, and the
collocation node family and count per element. Passed to
``FFModel.add_domain`` (and to ``add_input`` for an input's own node family).

Its enumerations double as the masks and sides used elsewhere:

- ``DOM``: ``ALL`` (the closed domain [a,b]), ``LB`` (a) and ``UB`` (b); masks
  combine by arithmetic, e.g. ``FFDom.ALL - FFDom.LB`` is (a,b];
- ``SIDE``: ``MINUS`` (tau^-, default) and ``PLUS`` (tau^+), for one-sided
  evaluations;
- ``TYPE``: ``LG``, ``LGR``, ``LGL``, ``CGL`` node families.

Example:

>>> from cronos import FFDom
>>> T = FFDom(0., 1., 10, FFDom.LGR, 4)  # 10 elements, 4 Radau nodes each
>>> Z = FFDom([0., .2, 1.], FFDom.LGL, 6)  # 2 elements from boundaries
)doc");

  // DOM is an enum.IntEnum, unlike pymcpp's py::enum_ elsewhere: masks combine
  // by subtraction (FFDom.ALL - FFDom.LB is (a,b]), which py::enum_ does not
  // support (it defines &, |, ^ only).
  py::native_enum<mc::FFDom::DOM>(pyFFDom, "DOM", "enum.IntEnum",
                                  "Domain masks for equations and outputs.")
      .value("ALL", mc::FFDom::ALL, "the closed domain [a,b]")
      .value("LB", mc::FFDom::LB, "the lower bound a")
      .value("UB", mc::FFDom::UB, "the upper bound b")
      .export_values()
      .finalize();
  py::enum_<mc::FFDom::SIDE>(pyFFDom, "SIDE", "Side of a point evaluation.")
      .value("MINUS", mc::FFDom::MINUS,
             "tau^-: the end of the element ending at tau (default)")
      .value("PLUS", mc::FFDom::PLUS,
             "tau^+: the start of the element starting at tau")
      .export_values();
  py::enum_<mc::FFDom::TYPE>(pyFFDom, "TYPE", "Collocation node family.")
      .value("LG", mc::FFDom::LG, "Legendre-Gauss")
      .value("LGR", mc::FFDom::LGR,
             "Legendre-Gauss-Radau, left endpoint included")
      .value("LGL", mc::FFDom::LGL,
             "Legendre-Gauss-Lobatto, both endpoints included")
      .value("CGL", mc::FFDom::CGL,
             "Chebyshev-Gauss-Lobatto, both endpoints included")
      .export_values();

  // --- FFDom Exceptions ---
  py::class_<mc::FFDom::Exceptions> pyFFDomExceptions(pyFFDom, "Exceptions",
                                                      R"doc(
Exception information for errors raised while defining a domain.

Instances carry an error code (ierr) and a human-readable description (what).
)doc");
  py::enum_<mc::FFDom::Exceptions::TYPE>(pyFFDomExceptions, "TYPE",
                                         "Error codes for FFDom exceptions")
      .value("BOUNDS", mc::FFDom::Exceptions::BOUNDS, "Invalid domain bounds")
      .value("ELEMENTS", mc::FFDom::Exceptions::ELEMENTS,
             "Invalid number of elements")
      .value("NODES", mc::FFDom::Exceptions::NODES,
             "Invalid number of collocation points")
      .value("UNDEF", mc::FFDom::Exceptions::UNDEF, "Undefined")
      .export_values();
  pyFFDomExceptions
      .def("ierr", &mc::FFDom::Exceptions::ierr, R"doc(
Return the error code as an integer (see FFDom.Exceptions.TYPE).
)doc")
      .def("what", &mc::FFDom::Exceptions::what, R"doc(
Return a string describing the error.
)doc");

  // --- FFDom Main Class ---
  pyFFDom.def(py::init<>(), "Empty domain.")
      .def(py::init<double const&, double const&, size_t,
                    mc::FFDom::TYPE const&, size_t>(),
           py::arg("lo"), py::arg("up"), py::arg("n_elem") = 1,
           py::arg_v("type", mc::FFDom::LG, "FFDom.LG"), py::arg("n_node") = 4,
           "Uniform finite elements on [lo, up] (one element by default).")
      .def(py::init<std::vector<double> const&, mc::FFDom::TYPE const&,
                    size_t>(),
           py::arg("elem_bnd"), py::arg_v("type", mc::FFDom::LG, "FFDom.LG"),
           py::arg("n_node") = 4,
           "Finite elements from their boundaries (size n_elem+1).")
      .def(py::init<double const&, std::vector<double> const&,
                    mc::FFDom::TYPE const&, size_t>(),
           py::arg("lo"), py::arg("elem_len"),
           py::arg_v("type", mc::FFDom::LG, "FFDom.LG"), py::arg("n_node") = 4,
           "Finite elements from the lower bound and the element lengths.")
      .def(py::init<mc::FFDom const&>(), py::arg("other"), "Copy.")
      .def_readonly("lo_dom", &mc::FFDom::lo_dom, "Lower bound.")
      .def_readonly("up_dom", &mc::FFDom::up_dom, "Upper bound.")
      .def_readonly("type", &mc::FFDom::type, "Collocation node family.")
      .def_readonly("n_elem", &mc::FFDom::n_elem, "Number of finite elements.")
      .def_readonly("n_node", &mc::FFDom::n_node,
                    "Number of collocation nodes per element.")
      .def_readonly("elem_bnd", &mc::FFDom::elem_bnd,
                    "Element boundaries, size n_elem+1.")
      .def_readonly("elem_len", &mc::FFDom::elem_len,
                    "Element lengths, size n_elem.")
      .def("uniform", &mc::FFDom::uniform,
           "True if all elements have the same length (to tolerance).")
      .def("elem_lo", &mc::FFDom::elem_lo, py::arg("iel"),
           "Lower bound of element ``iel``.")
      .def("elem_up", &mc::FFDom::elem_up, py::arg("iel"),
           "Upper bound of element ``iel``.")
      .def("elem_width", &mc::FFDom::elem_width, py::arg("iel"),
           "Length of element ``iel``.")
      .def("__repr__",
           [](mc::FFDom const& d)
           {
             static char const* T[] = {"LG", "LGR", "LGL", "CGL"};
             std::ostringstream os;
             os << "FFDom([" << d.lo_dom << ", " << d.up_dom
                << "], n_elem=" << d.n_elem << ", type=" << T[(int)d.type]
                << ", n_node=" << d.n_node << ")";
             return os.str();
           });
}
