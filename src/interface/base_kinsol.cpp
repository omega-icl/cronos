// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

// Python binding of BASE_KINSOL (base_kinsol.hpp): the corrector behind the consistent
// initial-condition solve -- its knobs, its report and its status.  A binder of its OWN because the
// types are shared: ODESLV's and DAESLV's INITIAL groups both hold a Knobs by value, so whichever
// solver is bound first would otherwise have to register them, and the second would find them
// already there.  Registered here once, both binders just reference them.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "base_kinsol.hpp"

namespace py = pybind11;

void
mc_base_kinsol(py::module_& m)
{
  py::class_<mc::BASE_KINSOL::Knobs> pyKnobs(m, "KINSOLKnobs", R"doc(
Knobs of the consistent initial-condition solve, held by value in a solver's
``INITIAL`` option group (``ODESLV.Options.INITIAL.KIN``, and the same on the DAE
side). The solve poses the INITIAL rows together with the dynamics residual at t0
over the states and their time derivatives, so these are Newton's knobs, not the
integrator's.

>>> o = ODESLV.Options()
>>> o.INITIAL.KIN.FTOL = 1e-12
>>> o.INITIAL.BYPASS_EXPLICIT = False   # force the solve on an explicit model
)doc");
  pyKnobs.def(py::init<>());
  pyKnobs.def_readwrite("FTOL", &mc::BASE_KINSOL::Knobs::FTOL,
                        "stopping tolerance on the scaled residual norm");
  pyKnobs.def_readwrite("STEPTOL", &mc::BASE_KINSOL::Knobs::STEPTOL,
                        "stopping tolerance on the scaled step -- what stops a slowly converging iteration");
  pyKnobs.def_readwrite("MAXIT", &mc::BASE_KINSOL::Knobs::MAXIT, "iteration cap");
  pyKnobs.def_readwrite("LINESEARCH", &mc::BASE_KINSOL::Knobs::LINESEARCH,
                        "use the line-search globalisation rather than plain Newton");
  pyKnobs.def_readwrite("MAXSETUP", &mc::BASE_KINSOL::Knobs::MAXSETUP,
                        "Jacobian refresh interval: 1 = every iteration");
  pyKnobs.def_readwrite("SPARSE", &mc::BASE_KINSOL::Knobs::SPARSE,
                        "use the sparse KLU linear solver rather than dense LU");
  pyKnobs.def_readwrite("SLOW_ITERS", &mc::BASE_KINSOL::Knobs::SLOW_ITERS,
                        "iteration count above which a CONVERGED solve is read as linearly convergent, "
                        "i.e. the rank-deficiency signature");

  py::enum_<mc::BASE_KINSOL::Status>(m, "KINSOLStatus", "Outcome of a consistent initial-condition solve.")
      .value("OK", mc::BASE_KINSOL::OK, "converged to FTOL")
      .value("AT_ROOT", mc::BASE_KINSOL::AT_ROOT,
             "the starting point already satisfied FTOL: nothing was solved")
      .value("RANK_DEFICIENT", mc::BASE_KINSOL::RANK_DEFICIENT,
             "converged linearly or stalled: the Jacobian is rank-deficient at this root")
      .value("FAILED", mc::BASE_KINSOL::FAILED, "KINSOL returned an error that is not the plateau signature")
      .export_values();

  py::class_<mc::BASE_KINSOL::Report> pyReport(m, "KINSOLReport", R"doc(
Outcome of one consistent initial-condition solve. ``RANK_DEFICIENT`` is a
diagnosis, not a failure: the iterate is accurate only to sqrt(FTOL) there, because
the residual behaves like the square of the error at a double root.
)doc");
  pyReport.def(py::init<>());
  pyReport.def_readonly("status", &mc::BASE_KINSOL::Report::status);
  pyReport.def_readonly("kin_flag", &mc::BASE_KINSOL::Report::kin_flag, "KINSOL's own return code");
  pyReport.def_readonly("iterations", &mc::BASE_KINSOL::Report::iterations);
  pyReport.def_readonly("fnorm", &mc::BASE_KINSOL::Report::fnorm, "final scaled residual norm");
  pyReport.def_readonly("message", &mc::BASE_KINSOL::Report::message);

}
