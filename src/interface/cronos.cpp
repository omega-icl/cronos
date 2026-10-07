// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

// Python module cronos: entry point.

#include <pybind11/pybind11.h>

namespace py = pybind11;

void mc_ffdom(py::module_&);
void mc_ffmodel(py::module_&);
void mc_odeslv(py::module_&);
void mc_ffodeslv(py::module_&);
void mc_ocfeslv(py::module_&);
void mc_ffocfe(py::module_&);

// What this binary was compiled with (the macros CMake sets on every source of the module).
static constexpr char const* kVersion = "5.0.0";
#if defined(CRONOS__WITH_KLU)
static constexpr bool kKLU = true;
#else
static constexpr bool kKLU = false;
#endif
#if defined(CRONOS__WITH_UMFPACK)
static constexpr bool kUMFPACK = true;
#else
static constexpr bool kUMFPACK = false;
#endif
#if defined(CRONOS__WITH_SPQR)
static constexpr bool kSPQR = true;
#else
static constexpr bool kSPQR = false;
#endif
#if defined(CRONOS__WITH_EIGEN)
static constexpr bool kEigen = true;
#else
static constexpr bool kEigen = false;
#endif
#if defined(ARMA_USE_SUPERLU)
static constexpr bool kSuperLU = true;
#else
static constexpr bool kSuperLU = false;
#endif

PYBIND11_MODULE(cronos, m)
{
  m.doc() = R"doc(
Python interface of CRONOS: models of differential-algebraic systems on the DAG of MC++ (``FFModel``), and their
solvers -- ``ODESLV`` (CVODES, with forward and adjoint sensitivities) and ``OCFESLV`` (orthogonal collocation on
finite elements) -- plus their embedding as DAG operations (``FFODESLV``, ``FFOCFESLV``, ``FFOCFERES``).

The DAG types (``FFGraph``, ``FFVar``, ``FFOp``) and the operators ``FFPartial``, ``FFIntegral``, ``FFEval`` come
from ``pymcpp``, which this module imports.
)doc";
  py::module_::import("pymcpp");  // registers FFBase/FFGraph/FFVar/FFOp: cronos
                                  // must share, not re-register them
  mc_ffdom(m);
  mc_ffmodel(m);
  mc_odeslv(m);    // after mc_ffmodel: ODESLV derives from FFModel
  mc_ffodeslv(m);  // after mc_odeslv: FFODESLV takes an ODESLV
  mc_ocfeslv(m);   // after mc_ffmodel: OCFESLV derives from FFModel
  mc_ffocfe(m);    // after mc_ocfeslv: FFOCFESLV / FFOCFERES take an OCFESLV
  
  m.attr("__version__") = kVersion;

  // build_info (2026-10-07): which build this binary is.  The licence follows from what is COMPILED IN, read from the
  // same macros the headers use, so it cannot disagree with the binary: UMFPACK and SPQR are GPL-2.0-or-later.
  m.def( "build_info", [](){
    py::dict b;
    b["KLU"]     = kKLU;
    b["UMFPACK"] = kUMFPACK;
    b["SPQR"]    = kSPQR;
    b["Eigen"]   = kEigen;
    b["SuperLU"] = kSuperLU;
    bool const gpl = kUMFPACK || kSPQR;
    py::dict d;
    d["version"]  = kVersion;
    d["license"]  = gpl ? "GPL-2.0-or-later" : "EPL-2.0";
    d["license_note"] = gpl
      ? "links GPL-2.0-or-later SuiteSparse components (UMFPACK and/or SPQR): this binary is distributed under the GPL"
      : "links no GPL component: this binary is distributed under the Eclipse Public License 2.0";
    d["backends"] = b;
    py::object pm = py::module_::import( "pymcpp" );
    py::object pv = py::none();                   // not a ?: -- that would convert the version string to None
    if( py::hasattr( pm, "__version__" ) ) pv = pm.attr( "__version__" );
    d["pymcpp"] = pv;
    return d;
  }, R"doc(
Describe this build of CRONOS: ``version``, ``license`` (``"EPL-2.0"``, or ``"GPL-2.0-or-later"`` when the GPL
SuiteSparse backends UMFPACK or SPQR are compiled in), ``license_note``, ``backends`` (which linear-algebra backends
are compiled in) and ``pymcpp`` (the version of the pymcpp module loaded alongside).
)doc" );
}
