# Building CRONOS with CMake

The CMake build mirrors MC++ v5's: a header-only interface target (`cronos_headers`, alias `cronos::headers`)
carrying the include paths, compile definitions and link libraries; the `cronos` Python extension module; its
`cronos.pyi` type stub; optional test drivers registered with CTest.

## Layout expected

```
<cronosroot>/CMakeLists.txt
<cronosroot>/src/*.hpp                      the 11 CRONOS headers
<cronosroot>/src/interface/*.cpp            the 7 binders (ffdom ffmodel odeslv ffode ocfeslv ffocfe cronos)
<cronosroot>/src/interface/gen_*_options.hpp  the 4 generated option headers the binders include
<cronosroot>/extern/pybind11                optional vendored pybind11 (else an installed package is used)
<cronosroot>/test/ODESLVS/*.cpp             ODESLV drivers
<cronosroot>/test/OCFESLV/*.cpp             OCFESLV drivers (+ helpers such as test_deriv_utils.hpp)
<cronosroot>/notebook/                      optional: tutorials, installed as-is
```

## Configure, build, install

```
cmake -S <cronosroot> -B build \
      -DMCPP_ROOT=<MC++ source tree, holding src/mc> \
      -DSUNDIALS_DIR=<sundials prefix>/lib/cmake/sundials \
      -DSUPERLU_ROOT=<SuperLU prefix>            # if not on the default search path
cmake --build build -j $(nproc)
cmake --install build --prefix <prefix>
```

`MCPP_ROOT`, `SUNDIALS_ROOT`, `SUITESPARSE_ROOT` and `SUPERLU_ROOT` are also read from the environment.

## Options

| option | default | effect |
|---|---|---|
| `MC_INTERVAL_LIBRARY` | `BOOST` | interval backend, as MC++: BOOST, PROFIL (`PROFIL_HOME`), FILIB (`FILIB_HOME`), NONVERIFIED |
| `ENABLE_HSL` | OFF | MC13/MC21/MC33 (`MC__USE_HSL`), as MC++ |
| `MC__USE_FADBAD` | OFF | FADBAD++ forward AD (reads `MCPP_ROOT/src/3rdparty/fadbad++`) |
| `CRONOS_WITH_KLU` | ON | sparse KLU in CVODES (`CRONOS__WITH_KLU`): SuiteSparse KLU + SUNDIALS built with KLU |
| `CRONOS_WITH_UMFPACK` | **OFF** | UMFPACK sparse LU in OCFESLV's setup (`CRONOS__WITH_UMFPACK`) -- **GPL-2.0-or-later**: a binary linking it is subject to the GPL; without it the setup uses a dense inverse |
| `CRONOS_WITH_EIGEN` | ON | Eigen sparse QR in OCFESLV (`CRONOS__WITH_EIGEN`) |
| `CRONOS_WITH_SUPERLU` | ON | SuperLU through `arma::spsolve` (`ARMA_USE_SUPERLU`), OCFESLV's default factorization |
| `CRONOS_WITH_SPQR` | **OFF** | SPQR factorization in OCFESLV (`CRONOS__WITH_SPQR`) -- **GPL-2.0-or-later**: a binary linking it is subject to the GPL |
| `CRONOS_ARMA_WRAPPER` | OFF | link Armadillo's runtime wrapper library instead of using Armadillo header-only with BLAS/LAPACK (and SuperLU) linked directly; the wrapper also pulls in whatever Armadillo was built with, e.g. ARPACK and MPI |
| `ENABLE_PYTHON` | ON | the `cronos` module |
| `CRONOS_STUBS` | ON | generate and install `cronos.pyi` (needs `pybind11-stubgen` and `pymcpp` importable; `PYMCPP_DIR` points at the directory holding `pymcpp` if it is not on the Python path) |
| `ENABLE_EXAMPLES` | OFF | build every driver in `test/ODESLVS` and `test/OCFESLV`, register each with CTest (labels `ODESLVS`, `OCFESLV`) |
| `CRONOS_CXX_STANDARD` | 17 | as the Make setup; 20 accepted |
| `CUSTOM_PYTHON_PATH` | | Python executable or virtual-environment root, as MC++ |

Each library is found through its own CMake package where it installs one (SUNDIALS >= 6, SuiteSparse >= 7,
Eigen) and otherwise by a library and header search, so older SuiteSparse releases and hand-installed SUNDIALS
work too.

## Python

* `cronos` uses `pymcpp`'s types (FFVar, FFGraph, FFPartial, ...) across the module boundary: **both modules must be
  built with the same pybind11 release** (v3.0.4 today: v3.1.0 has the `keep_alive<0,N>` regression on rejected
  overloads) and the same compiler; rebuild both after a pybind11 or header change.
* Installed into `${CMAKE_INSTALL_PREFIX}/lib` by default; for wheel builds with scikit-build-core, pass
  `-DCRONOS_INSTALL_DIR=.` and install only the `python_modules` component.

## Python checks: `make check`

```
cmake --build build --target check        # or: make -C build check
```
builds the `cronos` module if needed, then runs `check_ffmodel.py`, `check_odeslv.py`, `check_ffodeslv.py`,
`check_ocfeslv.py` and `check_ffocfe.py` in turn, stopping at the first failure.  They are taken from
`CRONOS_PYCHECK_DIR` (default `src/interface`) and find the modules through `CRONOS_PYPATH`, which the target sets
to the build directory followed by `PYMCPP_DIR`, so `pymcpp` must be importable from there (or from the Python
path).  Both modules must have been built from the same MC++ headers and pybind11 release.

## Tests

```
cmake -S <cronosroot> -B build -DENABLE_EXAMPLES=ON ...
cmake --build build -j $(nproc)
ctest --test-dir build -LE slow --output-on-failure          # the quick set
ctest --test-dir build -L OCFESLV -j4 --output-on-failure    # one directory (labels: ODESLV, OCFESLV)
ctest --test-dir build -L slow --output-on-failure           # the long drivers only
```
Each driver `X` is built from `X.cpp` in `test/ODESLV` (or `test/ODESLVS`) and `test/OCFESLV`, as the Make setup;
the OCFESLV drivers are compiled with `OCFE_OCFESLV_HEADER="ocfeslv.hpp"`.  Drivers taking 30 s or more on the
3 Oct 2026 sweep are labelled `slow` (`CRONOS_SLOW_TESTS`: PDE5, PDE5b, PSA7, PSA10, MBC3, MBC4, the MBC5/MBC6
family) with a timeout of `CRONOS_SLOW_TEST_TIMEOUT` (1800 s); the others time out after `CRONOS_TEST_TIMEOUT`
(300 s).

## MC++

MC++ v5's CMake installs headers but exports no package, so CRONOS reads them from `MCPP_ROOT/src/mc`.  The
CRONOS header set ships an updated `ocbase.hpp` (FFDom): it belongs in MC++ and must be pushed there.
