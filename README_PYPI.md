# cronos-mcpp: simulation and sensitivity analysis of ODE, DAE and PDE models in Python

`cronos-mcpp` is the Python distribution of [CRONOS](https://github.com/omega-icl/cronos), and imports as `cronos`.
It simulates models described by algebraic, ordinary differential, differential-algebraic and partial-differential
equations, computes the sensitivities of their outputs -- forward and adjoint -- and turns a simulation into a
differentiable function that can be composed with others, for optimisation and optimal control.

Models are declared on the factorable-function graphs of [pymcpp](https://pypi.org/project/pymcpp/), MC++'s Python
library, which supplies the symbolic representation and the automatic differentiation.  The solvers are SUNDIALS
CVODES, IDAS and KINSOL, and orthogonal collocation on finite elements, all compiled in C++.

## Installation

```
pip install cronos-mcpp
```

Wheels are built for Linux (x86_64), macOS (arm64) and Windows (x86_64), CPython 3.10 and later.  `pymcpp` and NumPy
are installed with it; `pymcpp` must be of the same release line, which the dependency enforces.

This is the **EPL build**, which links no GPL component.  A **GPL build**, with the UMFPACK and SPQR sparse backends,
is attached to each [GitHub release](https://github.com/omega-icl/cronos/releases) under the same name; giving pip the
release page selects it:

```
pip install cronos-mcpp -f https://github.com/omega-icl/cronos/releases/expanded_assets/<tag>
```

See [INSTALL.md](https://github.com/omega-icl/cronos/blob/5.x/INSTALL.md) for the differences between the two and
for building from source.

## What it does

* **Model declaration** (`FFModel`): states, inputs and constants on a pymcpp graph; one or several domains -- time
  and space -- discretised into finite elements; equations with their role (interior, initial, boundary, ...);
  point, distributed and integral outputs; transitions at given times; automatic order reduction of PDEs and index
  reduction of DAEs.
* **ODE integration** (`ODESLV`), with SUNDIALS CVODES.
* **DAE integration** (`DAESLV`), with SUNDIALS IDAS: implicit or singular mass matrices, consistent initialisation
  including the hidden constraints of high-index models, so only the free initial data are declared.
* **PDE models by the method of lines** (`DAESLV`): a model distributed in space is collocated in space and the
  resulting DAE integrated in time -- from the same declaration as `OCFESLV` takes.
* **PDE and DAE models by space-time collocation** (`OCFESLV`): orthogonal collocation on finite elements, solved as
  one system or by marching through the time elements.
* **Sensitivities**, forward and adjoint, of the outputs with respect to the inputs, set and read by name.
* **Embedding** (`FFODESLV`, `FFDAESLV`, `FFOCFESLV`): a solver as an operation of a pymcpp graph, to compose,
  evaluate in parallel and differentiate.

Results come back as NumPy arrays, every solver's options carry their documentation, and a `.pyi` stub serves
editors.

## An example

Heat conduction in a rod, `u_t = 0.1 u_zz` on `0 < z < 1`, held at a value `p` at `z = 0` and insulated at `z = 1`,
with the mean temperature at `t = 1` as output.  The declaration is written once; `DAESLV` solves it by the method of
lines and `OCFESLV` by space-time collocation:

```python
from pymcpp import FFGraph, FFPartial, FFIntegral
from cronos import FFDom, FFModel, DAESLV, OCFESLV

OpP, OpI = FFPartial(), FFIntegral()
Role = FFModel.EqnRole
io, ib, ii = FFModel.EqnOptions(Role.INTERIOR), FFModel.EqnOptions(Role.BOUNDARY), FFModel.EqnOptions(Role.INITIAL)
TI = FFDom.ALL - FFDom.LB                    # t > 0
ZI = FFDom.ALL - FFDom.LB - FFDom.UB         # 0 < z < 1

def rod(Solver):
    G = FFGraph()
    t, z, u, p = G.add_var("t"), G.add_var("z"), G.add_var("u"), G.add_var("p")
    S = Solver(G)
    S.add_domain(t, FFDom(0., 1., 2, FFDom.LGR, 3))    # in t: 2 elements
    S.add_domain(z, FFDom(0., 1., 3, FFDom.LGL, 5))    # in z: 3 elements of 5 nodes
    S.set_evolution_domain(t)
    S.add_state(u, [t, z], ref=0.5)
    S.add_input(p, ref=1.)
    S.add_equation(OpP(u, t) - OpP(0.1 * OpP(u, z), z), [t, z], [TI, ZI], io)
    S.add_equation(u - p,     [t, z], [TI, FFDom.LB], ib)
    S.add_equation(OpP(u, z), [t, z], [TI, FFDom.UB], ib)
    S.add_equation(u - p * (1. - z) * (1. - z), [t, z], [FFDom.LB, FFDom.ALL], ii)
    S.add_output(OpI(u, z), [t], point=[1.])          # the mean of u(1, z)
    return S, p

S, p = rod(DAESLV)                                     # rod(OCFESLV): space-time collocation
S.setup()
S.solve_fsens({p: 0.7})                                # or solve_asens: the adjoint
print("J =", S.val_function()[0], "  dJ/dp =", S.val_function_gradient()[0][0])
```

## Tutorials

Jupyter notebooks, each also a Python script, in the repository's
[notebook](https://github.com/omega-icl/cronos/tree/5.x/notebook) folder:

* [FFModel_tutorial](https://github.com/omega-icl/cronos/blob/5.x/notebook/FFModel_tutorial.ipynb) -- analysing a
  model before solving it: classification of the equations, index and order reduction, the model report;
* [ODESLV_tutorial](https://github.com/omega-icl/cronos/blob/5.x/notebook/ODESLV_tutorial.ipynb) -- a batch-reactor
  network: ODE simulation, forward and adjoint sensitivities;
* [DAESLV_tutorial](https://github.com/omega-icl/cronos/blob/5.x/notebook/DAESLV_tutorial.ipynb) -- a vented batch
  reactor, an index-2 DAE: consistent initialisation, sensitivities, parallel scenario evaluation;
* [OCFESLV_tutorial](https://github.com/omega-icl/cronos/blob/5.x/notebook/OCFESLV_tutorial.ipynb) -- nonlinear heat
  conduction in a rod: both solve modes, sensitivities, the embedding in a graph;
* [MOL_tutorial](https://github.com/omega-icl/cronos/blob/5.x/notebook/MOL_tutorial.ipynb) -- the same rod by the
  method of lines, compared with space-time collocation;
* [PSA_tutorial](https://github.com/omega-icl/cronos/blob/5.x/notebook/PSA_tutorial.ipynb) -- breakthrough in a
  non-isothermal adsorption column, by either solver.

## Documentation and source

CRONOS is also a header-only C++ library; the source, its documentation, the release notes
([CHANGES.md](https://github.com/omega-icl/cronos/blob/5.x/CHANGES.md)) and the issue tracker are on
[GitHub](https://github.com/omega-icl/cronos).

## Licence

The EPL build is distributed under the Eclipse Public License 2.0; the GPL build under GPL-2.0-or-later, as a
Secondary License.  The third-party components are listed in THIRD_PARTY_NOTICES.txt.  CRONOS is developed by
Benoit Chachuat and the OMEGA Research Group at Imperial College London.
