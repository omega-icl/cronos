# Installing CRONOS

## What CRONOS is

CRONOS is a header-only C++ library for the simulation and sensitivity analysis of dynamic models -- ODEs, DAEs
and PDE/DAE systems -- and for embedding those simulations in optimisation and optimal control.  It lives in the
`mc::` namespace and is built on [MC++](https://github.com/omega-icl/mcpp), whose directed acyclic graphs (DAGs)
of factorable functions provide the model representation and automatic differentiation.

* **`FFModel`** declares a model on an MC++ DAG: states, inputs and constants, their domains (time and space,
  discretised by finite elements), the equations with their role (interior, initial, boundary, ...), and the
  outputs.
* **`ODESLV`** (`odeslvs_cvodes.hpp`) integrates ODE/DAE models with SUNDIALS CVODES, with forward and adjoint
  sensitivities.
* **`OCFESLV`** (`ocfeslv.hpp`) solves PDE/DAE models by orthogonal collocation on finite elements, either as one
  monolithic system or by marching through the time elements, with forward and adjoint sensitivities and state
  transitions (jumps); the monolithic time discretisation is causal, so both modes solve the same discretisation.
* **`FFODESLV`, `FFOCFESLV`, `FFOCFERES`** embed a solver -- or the collocation residuals -- in an MC++ DAG as an
  external operation, so a simulation becomes a DAG function that can be composed, evaluated, sampled in parallel
  and differentiated.

A Python module, `cronos`, exposes all of this on top of PyMC++, MC++'s Python library (module `pymcpp`).

## Requirements

| dependency | version | role |
|---|---|---|
| C++ compiler | C++17 | developed and tested with GCC |
| CMake | >= 3.18 | |
| MC++ | v5 source tree | headers read from `<MC++>/src/mc` |
| Boost | | intervals (with MC++'s default interval backend) |
| BLAS / LAPACK | | |
| Armadillo | | dense and sparse linear algebra |
| SuiteSparse | 7.x recommended (older found by library search) | KLU (LGPL-2.1+); UMFPACK and SPQR optional (GPL-2.0-or-later, see below) |
| SuperLU | 5.2.x | OCFESLV's default sparse factorisation, through Armadillo |
| Eigen | >= 3.3 | optional sparse QR |
| SUNDIALS | **>= 7**, built with KLU (`-DENABLE_KLU=ON`) | CVODES for ODESLV |
| Python | 3.10+ (development headers) | for the `cronos` module |
| pybind11 | **>= 3.0.3, not 3.1.0** (3.0.4 recommended) | Python bindings |
| pybind11-stubgen | optional | generates the `cronos.pyi` type stub |
| HSL MC13/MC21/MC33 | optional | as for MC++ |

Two compatibility notes from experience:
* **pybind11 3.1.0** has a `keep_alive` regression on rejected overloads that breaks the bindings; the build refuses
  it.
* **Armadillo 12 with SuperLU 6** is not a supported combination; use SuperLU 5.2.x.

## 1. Build MC++ and `pymcpp` first

CRONOS reads MC++'s headers, and the `cronos` module passes MC++ objects (`FFVar`, `FFGraph`, ...) to and from
`pymcpp`.  pybind11 shares types between two modules only when both were built against the same internals, so:

* build `pymcpp` with **the same pybind11 release** (3.0.4) as `cronos`;
* build both from **the same MC++ headers**, and rebuild both after any MC++ header change.  A mismatch typically
  shows as `free(): invalid pointer` or a segmentation fault as soon as the two modules exchange objects.

The CRONOS header set may ship an updated MC++ header (currently `ocbase.hpp`): copy it into `<MC++>/src/mc` before
building either.

## 2. Configure

```
cmake -S <cronos> -B build \
      -DMCPP_ROOT=<MC++ source tree> \
      -DPYMCPP_DIR=<directory holding the pymcpp module>
```

Most dependencies are found without help when installed in standard places.  Otherwise:

| variable (cache or environment) | points at |
|---|---|
| `MCPP_ROOT` | the MC++ source tree (required) |
| `SUNDIALS_ROOT` | the SUNDIALS prefix -- not needed if SUNDIALS is in `/usr`, `/usr/local`, `/opt/sundials*` or on `CMAKE_PREFIX_PATH` |
| `SUITESPARSE_ROOT` | the SuiteSparse prefix |
| `SUPERLU_ROOT` | the SuperLU prefix |
| `PYMCPP_DIR` | the directory holding `pymcpp`, if it is not on the Python path |
| `CUSTOM_PYTHON_PATH` | a Python executable or virtual environment, as for MC++ |

pybind11 is taken from `extern/pybind11` if present, else from an installed pybind11 package, else v3.0.4 is fetched
from GitHub (`-DCRONOS_FETCH_PYBIND11=OFF` to forbid it).  An unacceptable release is refused with the fix.

Main options (all listed in `BUILD_CMAKE.md`):

| option | default | |
|---|---|---|
| `CRONOS_WITH_KLU`, `CRONOS_WITH_EIGEN`, `CRONOS_WITH_SUPERLU` | ON | linear-solver backends |
| `CRONOS_WITH_UMFPACK` | **OFF** | sparse LU for the setup's trace elimination (else dense) -- **GPL-2.0-or-later**: a binary linking it is subject to the GPL |
| `CRONOS_WITH_SPQR` | **OFF** | SPQR factorisation, and the W-test's null bases above the dense cap -- **GPL-2.0-or-later**: a binary linking it is subject to the GPL |
| `MC_INTERVAL_LIBRARY` | BOOST | as MC++ (BOOST, PROFIL, FILIB, NONVERIFIED) |
| `ENABLE_HSL`, `MC__USE_FADBAD` | OFF | as MC++ |
| `ENABLE_PYTHON` | ON | the `cronos` module and its stub |
| `ENABLE_EXAMPLES` | OFF | the test drivers, registered with CTest |

## 3. Build and check

```
cmake --build build -j $(nproc)
cmake --build build --target check        # the five Python check scripts
```

`check` runs `check_ffmodel.py`, `check_odeslv.py`, `check_ffodeslv.py`, `check_ocfeslv.py` and `check_ffocfe.py`
from `src/interface` against the freshly built module, stopping at the first failure.  Each prints its checks and
a final `... passed, 0 failed -- ALL PASS`.

## 4. Test drivers (optional)

```
cmake -S <cronos> -B build -DENABLE_EXAMPLES=ON ...
cmake --build build -j $(nproc)
ctest --test-dir build -LE slow --output-on-failure       # the quick set
ctest --test-dir build -L slow --output-on-failure        # the long drivers (minutes each)
```

Each driver in `test/ODESLV` and `test/OCFESLV` is one executable and one CTest test, labelled by directory; the
drivers taking 30 s or more are also labelled `slow`.

## 5. Install

```
cmake --install build --prefix <prefix>
```

installs the headers in `<prefix>/include`, the `cronos` module and `cronos.pyi` in `<prefix>/lib`
(`CRONOS_INSTALL_DIR`), and the notebooks in `<prefix>/notebook`.  Use `cronos` from Python with `<prefix>/lib` and
`pymcpp`'s directory on the Python path:

```python
import pymcpp, cronos
```

The tutorials (`ODESLV_tutorial`, `OCFESLV_tutorial`) are the best starting point.

## Troubleshooting

| symptom | cause |
|---|---|
| `free(): invalid pointer` or a crash when `cronos` and `pymcpp` exchange objects | the two modules were built from different MC++ headers or pybind11 releases: rebuild both from clean |
| `SUNDIALS x.y found ..., but CRONOS requires SUNDIALS >= 7` | install SUNDIALS 7 and set `SUNDIALS_ROOT` |
| `... no target for component 'sunlinsolklu'` | SUNDIALS built without KLU: rebuild it with `-DENABLE_KLU=ON`, or `-DCRONOS_WITH_KLU=OFF` |
| `extern/pybind11 is pybind11 3.1.0` | check out v3.0.4 there, or remove `extern/pybind11` |
| `pybind11 changed from X to Y since this build directory was configured` | CMake cannot see a pybind11 change under the same path and would mix objects of two releases: start from a clean build directory, and rebuild `pymcpp` with the same release |
| the `cronos.pyi` stub is not generated | `pybind11-stubgen` or `pymcpp` not importable by the build's Python: `pip install pybind11-stubgen`, set `PYMCPP_DIR` |
| a stale option or docstring in Python | the `gen_*_options.hpp` headers are generated from the C++ headers (`gen_options.py`); regenerate with `cmake --build <dir> --target gen-options` (never edit; keep a single copy, next to the binders); `make check` verifies they are up to date |
| oversubscribed cores in a parallel `FFGraph::veval` over an embedded solver | set the solver's `options.MAXTHREAD = 1`: `veval` already runs the solves in parallel |

## Installing from PyPI

`pip install cronos-mcpp` installs the **EPL build**: wheels for Linux (x86_64) and macOS (arm64), CPython 3.10 and
later; the module imports as `cronos`, and its dependencies (`pymcpp`, NumPy) come from PyPI.  It links no GPL
component: OCFESLV eliminates the trace multipliers at setup with a dense inverse (fine up to a few thousand; cubic
beyond), and the SPQR-based options (`SOLVE_SPQR`, the W-test's null basis above its dense cap) are not available.

### The GPL build

With UMFPACK (sparse elimination at setup) and SPQR, under GPL-2.0-or-later.  Its wheels are attached to each GitHub
release under the same name, with a build tag; giving pip (or uv) the release page makes it prefer them:

```
pip install cronos-mcpp -f https://github.com/omega-icl/cronos/releases/expanded_assets/v5.0.0
uv pip install cronos-mcpp --find-links https://github.com/omega-icl/cronos/releases/expanded_assets/v5.0.0
```

The page must be given at EVERY install or upgrade: without it, pip finds only the PyPI wheels and installs the EPL
build.  To make that permanent, use one of: a line `-f https://github.com/omega-icl/cronos/releases/expanded_assets/v5.0.0`
in `requirements.txt`; `find-links` in pip's configuration (`pip config set global.find-links <page>`); or
`find-links = ["<page>"]` under `[tool.uv]` in a uv project.  A package that depends on `cronos-mcpp` is satisfied by
either build, so installing it keeps the GPL build.

To build it from source instead (the native dependencies above installed, and `MCPP_ROOT` pointing to MC++'s sources):

```
pip install cronos-mcpp --no-binary cronos-mcpp \
    -C cmake.define.CRONOS_WITH_UMFPACK=ON -C cmake.define.CRONOS_WITH_SPQR=ON
```

Check which build is installed:

```
python -c "import cronos; print(cronos.build_info())"     # 'license': 'EPL-2.0' or 'GPL-2.0-or-later'
```

## Licence

CRONOS is published under the Eclipse Public License 2.0, with the GNU General Public License, version 2 or later, as
a Secondary License (see [LICENSE](LICENSE)).  Two optional SuiteSparse backends are GPL-2.0-or-later: UMFPACK
(`CRONOS_WITH_UMFPACK`: sparse elimination of the trace multipliers at setup; without it a dense inverse is used,
slower on large interface systems) and SPQR (`CRONOS_WITH_SPQR`: a sparse QR factorisation, and the null bases of
the interface plan's W-test above its dense size cap).  Both are OFF by default, so a default build links no GPL
component; a binary that links either is subject to the GPL.  The `cronos-mcpp` wheels on PyPI are the EPL build;
the GPL build, with both backends, is published with each GitHub release (see [INSTALL.md](INSTALL.md));
`cronos.build_info()` reports which build is installed.
