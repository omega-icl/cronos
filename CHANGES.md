# Changes

## 5.1.0 (2026-10-27)

The headline is a second integrator. **DAESLV** integrates the model's undivided form with SUNDIALS IDAS,
beside ODESLV on CVODES, with the same interface, the same sensitivities, and a shared base. Around it: the
consistent-initialisation solve both integrators now take, index reduction by dummy derivatives as an
alternative to substitution, and a setup and solve path that scale to meshes of thousands of states.

### DAESLV: DAE integration on IDAS

* **`DAESLV_IDAS`** (`daeslv_idas.hpp`, with `daeslv_base.hpp` and `base_idas.hpp`): the model's UNDIVIDED form
  -- a mass matrix and a residual, so an implicit or singular mass matrix needs no division -- with the stage
  structure, transitions and recording of ODESLV, and dense or sparse (KLU) linear algebra.
* **Consistent initialisation** by a shared KINSOL corrector (`base_kinsol.hpp`): one parameterised system per
  role (initial point, stage restart, stop), with the high-index hidden constraints materialised as rows at
  the initial point (`REDUCE.HIDDEN_IC`, default on), so a model declares only its FREE initial data.  Both
  integrators take it; the explicit bypass sits behind `INITIAL.BYPASS_EXPLICIT`.
* **Forward sensitivities** through IDAS, seeded from the solved consistent-initialisation system, including
  across stops and transitions; quadrature and output sensitivities.
* **Adjoint sensitivities** (`solve_asens`): the terminal condition posed on the differential block and the
  seed -- at the final time and at every jump -- solved as one linear extended system by structure, with no
  `IDACalcICB` Newton iteration; the backward problems of SEVERAL outputs registered together and advanced in
  ONE sweep.  The backward residual carries the `(dF_y'/dt)^T lambda` term, so a state-dependent mass matrix
  is handled; a point output on an algebraic state is seeded by projection through the residual at the final
  time (Cao, Li & Petzold), with its direct term in the gradient.  `DAESLV_asensmulti` gates the adjoint
  against the forward on several outputs.
* **Transitions** as rows of the stop system; **outputs** with integral captures as quadratures and point
  captures as function terms; **recording** opt-in, per block, sub-sampled on the adjoint sweep, with
  `RESRECORD` meaning the same count in both solvers.
* **`SOLVE.VERBOSE`** honoured, a stage summary at every boundary; **`SOLVE.LINSOL`** DENSE or SPARSE.
* **Embedding**: `FFDAESLV` (`ffde.hpp`, formerly `ffode.hpp`, which now defines both operations as one
  template instantiated per solver) -- the DAE solver as an external operation of an MC++ DAG, composable,
  differentiable forward and reverse, and evaluable in parallel (`DAESLV_threads` gates the identical-variable
  sweep with two workers).  Requires MC++ with the `FFSubgraph::save_val()` slot guard.
* **Python**: `cronos.DAESLV_IDAS`, `cronos.FFDAESLV`, the KINSOL knobs, `check_daeslv.py` and
  `check_ffdaeslv.py`; the **`DAESLV_tutorial`** notebook and script -- a vented batch reactor, index 2, with a
  closed form to check against.
* **Tests**: 26 drivers in `test/DAESLV`, from the first integration to the adjoint on several outputs,
  threads, scaling to 9,600 states, and the two reductions side by side (`DAESLV_reduction`).

### Index reduction

* **Dummy derivatives** (`REDUCE.DUMMY_DERIV`, default off): a high-index block reduced by the method of
  Mattsson and Soderlind rather than by substitution -- every hidden level formed first, then as many
  derivatives chosen to be dummies as there are levels, under the nonsingularity test on their coefficient
  matrix at the reference values, preferring a state no output reads.  The constraint stays a row and holds
  to solver tolerance: 3e-13 on the index-3 pendulum against substitution's 2e-9.  Index 2 and 3; forward
  solves, forward and adjoint sensitivities, and transitions (the algebraic states re-solve at the stop) on
  both.  The static pivot is diagnosed: a solve
  stops where it degenerates, naming the time, the values and the remedies, rather than continuing to a
  wrong answer (`DAESLV_dummy`, `DAESLV_drift`).
* Substitution remains the default.  On a LINEAR constraint it holds the constraint to round-off (BDF
  conserves a linear invariant); dummy derivatives hold it to Newton's tolerance.  The drift substitution
  shows on NONLINEAR constraints is proportional to the integration tolerance.

### ODESLV and the shared base

* **`BASE_DE`** (formerly `BASE_DYNAMIC`; `FFBaseDE` formerly `FFBaseODESLV`): the stage machinery, results
  and recording, `fdiff`, `var_parameter`, the controls hook and the consistent-initialisation plumbing hoisted
  from ODESLV so both integrators share one body.
* **ODESLVS merged into ODESLV**: `odeslvs_base.hpp` and `odeslvs_cvodes.hpp` are gone; `ODESLV_CVODES` carries
  the sensitivities.  Options restructured into `SOLVE` / `OUTPUT`, mirroring OCFESLV and DAESLV.
* Fixes: the IC solve read its parameters one slot late; the forward sensitivity IC now comes from the solved
  system, not the divided form; the IC seed was the origin, not the declared reference; the base `record()`
  body and `using FFModel::fdiff` lost in the merge restored.

### Performance

Setup and solve were measured on a reaction-diffusion model by the method of lines, from 118 to 9,598 states
(`DAESLV_scaling`), and every O(n^2) step found was removed:

* **Setup**: extraction was n_eqn x n_x symbolic differentiations, now one subgraph walk per row; a time-only
  block's classification skips the principal symbol; every setup-time Jacobian -- the IC systems, dF/dy, the
  stop system's sensitivity seeds -- is built row by row (`_sfad_by_row`), and dF/dy' is assembled from the
  mass coefficients.  Setup at 2,400 states: 11 s to 1.3 s.
* **Solve**: `explicit_initial_conditions()` ran two dense n x n Jacobians at EVERY solve, and ODESLV's
  Jacobian setup one O(n^2) SFAD; both are structural and cached now.  IDAS's residual callback evaluated the
  residual AND every Jacobian entry; it evaluates the residual.  ODE at 9,600 states, sparse: 9.4 s to
  0.59 s.
* The DAE's remaining cost against the ODE is the index-1 formulation's own: more rows, more steps, more
  Newton work per step, in IDAS proper.

### FFModel and OCFESLV

* `add_equation` takes an `EqnRole` directly (C++ and Python); the classification provisions inputs and
  samples them with the states; one `[classify]` line per block at `DISPLAY_LEVEL >= 1`.
* **`CRONOS_DISPLAY_LEVEL`** is a DEFAULT, not a floor: it applies where a driver has not set `DISPLAY_LEVEL`,
  and a value the driver sets, `0` included, stands.
* Interface plan: a per-node receiver rule inferred by FFInv with `EqnOptions::defines` as fallback;
  `INTERFACE.DROP_POLICY` defaults to `DROP_OFF`.  Classification keeps the last provider of a derivative in
  the symbol; an unevaluable principal symbol is contained and the homotopy ramp saves the trace block; the
  W-test's verdict is independent of the audit and of stale caches.  Moving-mesh driver defaults revised.
* Python: `add_homotopy`, `reset_homotopy`, `set_homotopy_cap`, `continuation_report`; `eval_solution` takes a
  side; `setup_status_str`.

### Build and packaging

* Version 5.1.0; `pyproject.toml` fixes for the Windows wheel; the macOS wheel import fixed and the
  sensitivity vector arrays freed.
* `cmake --build --target check` runs the Python check scripts, `gen_options.py --check`,
  `gen_environment.py --check` and Doxygen with no warning allowed; the option headers are generated from
  the C++ structs and the environment-variable table from a registry.
* `.clang-format` identical to MC++; the DAESLV stack documented in `BUILD_CMAKE.md` and the headers.

### Known gaps

* The dummy pivot is static.  A trajectory through its singularity -- a rotating system, the pendulum's full
  swing -- stops with a diagnostic naming the time, the values and the remedies; re-selection per stage is
  the designed next step.
* Dummy derivatives decline, loudly and falling back to substitution, for a dummy that would be the block's
  only differential state, and for a level holding the derivative of a state with no `d/dt` row that is
  differentiated again (several dummies per level).

## 5.0.0 (2026-10-09)

`FFModel`, `ODESLV` on CVODES, `OCFESLV` with the interface-plan refactor, the DAG embeddings `FFODESLV`,
`FFOCFESLV` and `FFOCFERES`, and the `cronos` Python module on PyMC++.
