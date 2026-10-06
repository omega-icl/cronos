# CRONOS: Simulation and Sensitivity Analysis of Dynamic Models on Factorable-Function DAGs

CRONOS is a header-only C++ library for the simulation and sensitivity analysis of dynamic models -- systems of
ordinary differential equations (ODEs), differential-algebraic equations (DAEs) and partial differential equations
(PDE/DAE systems) -- and for embedding those simulations in optimisation and optimal control.  It is written in
C++ for execution speed, with Python binders in the module `cronos` through [pybind11](https://pybind11.readthedocs.io/).

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

The tutorials in `notebook/` (`ODESLV_tutorial`, `OCFESLV_tutorial`) walk through a batch-reactor network and a
heated rod: model declaration, both solve modes, sensitivities, and the embedding in a DAG.

## Setting up CRONOS

CRONOS depends on MC++ (v5), Boost, BLAS/LAPACK, Armadillo, SuiteSparse, SuperLU, Eigen and SUNDIALS (>= 7),
and on Python with pybind11 for the `cronos` module.  Refer to [INSTALL.md](INSTALL.md) for the installation
guide and to [BUILD_CMAKE.md](BUILD_CMAKE.md) for the complete list of build options.

```
cmake -S . -B build -DMCPP_ROOT=<MC++ source tree> -DPYMCPP_DIR=<directory holding pymcpp>
cmake --build build -j $(nproc)
cmake --build build --target check        # the Python check scripts
cmake --install build --prefix <prefix>
```

## Layout

```
src/                 the CRONOS headers (header-only library)
src/interface/       the Python binders, the generated option headers, the check scripts
test/ODESLV/         ODESLV test drivers
test/OCFESLV/        OCFESLV test drivers (the corpus the solver is validated against)
notebook/            tutorials (Python scripts and Jupyter notebooks)
```

## Contacts

* Repo owner: [Benoit C. Chachuat](https://profiles.imperial.ac.uk/b.chachuat)
* OMEGA Research Group, Imperial College London

## Licence

CRONOS is published under the Eclipse Public License.  The optional SPQR backend (SuiteSparse) is GPL-2 licensed:
binaries that link it are subject to the GPL; it is off by default.
