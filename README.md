# CRONOS: Simulation and Sensitivity Analysis of Dynamic Models on Factorable-Function DAGs

CRONOS is a header-only C++ library for the simulation and sensitivity analysis of mathematical models described
by systems of algebraic, differential and/or partial-differential equations, and for embedding those simulations
in optimisation and optimal control.  It is written in C++ for execution speed, with Python binders in the module
`cronos` through [pybind11](https://pybind11.readthedocs.io/).

CRONOS extends [MC++](https://github.com/omega-icl/mcpp): models are declared on MC++'s directed acyclic graphs
(DAGs) of factorable functions, which provide the symbolic representation and the automatic differentiation; and a
CRONOS simulation can itself become an *external operation* of an MC++ DAG, so that the solution of a dynamic
system is a factorable function that MC++ -- and the libraries built on it, such as
[CANON](https://github.com/omega-icl/canon) for numerical optimisation and [MAGNUS](https://github.com/omega-icl/magnus)
for model development and analysis -- can compose, evaluate and differentiate.

## Capabilities

* **Model declaration** (`FFModel`): states, inputs and constants on an MC++ DAG; one or several domains -- time
  and space -- each discretised into finite elements with a collocation node family; equations with their role
  (interior, initial, boundary, interface, link, ...) and the region they hold on; point, distributed and integral
  outputs; state transitions (jumps) at given times; automatic index reduction of DAEs and order reduction of PDEs.
* **ODE/DAE integration** (`ODESLV`, header `odeslvs_cvodes.hpp`): with SUNDIALS CVODES, with forward and adjoint
  sensitivities of the outputs with respect to the model's parameters and controls.
* **Orthogonal collocation on finite elements** (`OCFESLV`, header `ocfeslv.hpp`): PDE/DAE models on tensor
  domains, solved either as one monolithic system or by marching through the time elements -- both solve the
  same discretisation, since the monolithic time discretisation is causal -- with strong, weak (penalty) or
  trace imposition of the interface conditions, forward and adjoint sensitivities, and state transitions.
* **Embedding in DAGs** (`FFODESLV`, `FFOCFESLV`, `FFOCFERES`): a solver, or the collocation residuals, as an
  external operation of an MC++ DAG: the simulation becomes a function that can be composed with other DAG
  operations, evaluated and sampled in parallel, and differentiated -- forward and reverse -- through the
  solver's own sensitivities.
* **Python** (`cronos`): all of the above from Python, on top of PyMC++ (MC++'s Python library, module `pymcpp`),
  with NumPy arrays for results and a `.pyi` type stub for editors.

The tutorials in `notebook/`, each a Jupyter notebook and a Python script:

* `FFModel_tutorial` -- analysing a model before solving it: classification of the equations, index and order
  reduction, initial data for high-index DAEs, and the model report;
* `ODESLV_tutorial` -- a serial batch-reactor network: ODE simulation, and forward and adjoint sensitivities;
* `OCFESLV_tutorial` -- nonlinear heat conduction in a rod with a heater schedule: both solve modes, sensitivities,
  and the embedding in a DAG.

## Configuration

CRONOS is configured through its options -- `FFModel::Options` for the model, `OCFESLV::Options` and the ODESLV
options for the solvers -- available in Python as `options` on each object, with their documentation.  A few
`CRONOS_*` environment variables remain: diagnostics that instrument a whole test sweep without editing any driver,
a thread cap, and three test switches.  None is needed to run a model; they are listed in
[docs/ENVIRONMENT.md](docs/ENVIRONMENT.md).

## Installing from PyPI

```
pip install cronos-mcpp                    # the EPL build; imports as `cronos`
```

The GPL build, with the faster SuiteSparse backends UMFPACK and SPQR for large models, is a one-line install from the
release page -- see [INSTALL.md](INSTALL.md#the-gpl-build).

## Setting up CRONOS

CRONOS depends on MC++ (v5), Boost, BLAS/LAPACK, Armadillo, SuiteSparse, SuperLU, Eigen and SUNDIALS (>= 7),
and on Python with pybind11 for the `cronos` module.  Refer to [INSTALL.md](INSTALL.md) for the installation
guide and to [BUILD_CMAKE.md](BUILD_CMAKE.md) for the complete list of build options.

```
cmake -S . -B build -DMCPP_ROOT=<MC++ source tree> -DPYMCPP_DIR=<directory holding pymcpp>
cmake --build build -j $(nproc)
cmake --build build --target check        # the generated-file checks and the Python check scripts
cmake --install build --prefix <prefix>
```

The API documentation is generated with [Doxygen](https://www.doxygen.nl): `cd docs && doxygen CRONOS.dox`, then
open `docs/html/index.html`.

## Layout

```
src/                 the CRONOS headers (header-only library)
src/interface/       the Python binders, the generated option headers, the check scripts
test/ODESLV/         ODESLV test drivers
test/OCFESLV/        OCFESLV test drivers (the corpus the solver is validated against)
notebook/            tutorials (Python scripts and Jupyter notebooks)
docs/                the Doxygen configuration (CRONOS.dox) and main page (CRONOS.txt); ENVIRONMENT.md
```

## Contacts

* Repo owner: [Benoit C. Chachuat](https://profiles.imperial.ac.uk/b.chachuat)
* OMEGA Research Group, Imperial College London

## Licence

CRONOS is published under the Eclipse Public License 2.0, with the GNU General Public License, version 2 or later, as
a Secondary License (see [LICENSE](LICENSE)).  Two optional SuiteSparse backends are GPL-2.0-or-later: UMFPACK
(`CRONOS_WITH_UMFPACK`: sparse elimination of the trace multipliers at setup; without it a dense inverse is used,
slower on large interface systems) and SPQR (`CRONOS_WITH_SPQR`: a sparse QR factorisation, and the null bases of
the interface plan's W-test above its dense size cap).  Both are OFF by default, so a default build links no GPL
component; a binary that links either is subject to the GPL.  The `cronos-mcpp` wheels on PyPI are the EPL build;
the GPL build, with both backends, is published with each GitHub release (see [INSTALL.md](INSTALL.md));
`cronos.build_info()` reports which build is installed.
