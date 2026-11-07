# Changes

## Unreleased

* **DAESLV: the consistent initialisation converges to the integrator's tolerance.** KINSOL's point met the
  residual tolerance but pinned an ill-conditioned algebraic unknown only to ~1e-12, which IDAS error-tests at
  `ATOL`; a polish pass with KINSOL's scaled-step criterion in the error-weighted norm `1/(RTOL|y|+ATOL)` now ends
  every IC solve, at t0 and at each stop (`BASE_KINSOL::Knobs::RTOL/ATOL/WSTEPTOL`). This is what IDACalcIC did.
* **MoL on DAESLV (6c): time-varying inputs.** Each spatial slot of a declared input over `t` becomes an input over
  `t` of the semi-discrete DAE with the declared collocation; declared DOFs and working parameters are matched as
  a bijection, so by-name values and gradient rows read as the declared DOFs. `ic_sens` now differentiates through
  the placeholders binding a time-varying input at `t0` (an INITIAL row `x(0) = u(t0)` had no sensitivity to `u`,
  in ODESLV and DAESLV alike), and such a row no longer counts as explicit for the bypass.
* **MoL on DAESLV (6c): transitions.** A declared transition is applied at its state's differential nodes, each
  component composed pointwise there (`OCFEMOL::transition_expressions`, the INITIAL rows' form: the integrator
  peels these rows, so no vectorised operation); the integrator's stop machinery re-solves the algebraic nodes and
  carries the sensitivities across. `FFMOLINI`/`FFGradMOLINI`, the vectorised Ini stream, retired for that reason.
* **DAESLV adjoint: one backward problem for every adjoined output.** The outputs' multipliers are stacked into
  one IDA instance with one step sequence and one block-diagonal matrix (KLU whenever there is more than one
  block or the forward factorises sparse), so the forward rows and Jacobian are evaluated once per backward
  residual rather than once per output. The backward Jacobian is analytic (it was built by difference quotients),
  and the quadrature sensitivity integrands read the memoised evaluation. The PSA column at 480 states, one
  parameter, two outputs: `solve_fsens` 4.7 s -> 2.6 s, `solve_asens` 6.1 s -> 3.9 s, against 0.6 s for the solve.
* **MoL: the consistent initialisation starts from the declared profile at every node.** An INITIAL row at an
  algebraic node (a Dirichlet node, or by default a Neumann node) is set aside, the node being held by its
  boundary row; `OCFEMOL::init_seed` skipped those rows, so the node started from its reference value, and a
  boundary row with several roots could be met on the wrong one. On the heated rod with `D(u) = D0 (1 + beta u)`
  the flux row was met at `1 + beta u < 0`, a negative diffusivity: `u(2,1) = -1.02` against OCFESLV's `2.79`.
  The set-aside rows now seed their nodes too; DAESLV reports, at `DISPLAY_LEVEL >= 1`, how far the consistent
  initialisation moved such a node from its declared value (`initial_consistency_defect`, unread before).
* **Notebooks: `MOL_tutorial` and `PSA_tutorial`.** The OCFESLV tutorial's heated rod on DAESLV by the method of
  lines, one declaration for both solvers: outputs, profiles, OCFESLV's spectral convergence in its time order,
  gradients by name in the shared layout. A non-isothermal two-component adsorption column, by MoL or by space-time
  collocation behind one switch: breakthrough, thermal front, sensitivities to the affinity, a cross-check.
* **Outputs that are NONLINEAR functions of horizon integrals and point evaluations**, on ODESLV, DAESLV and
  under MoL (refused before: one record or an affine combination). Before the model sets up, each distinct
  reduction in `t` such an output reads is lifted into a hidden STATE -- an integral into the quadrature state
  `q' = g`, `q(0) = 0`; an interior point value into a frozen state with the transition `s(tau+) = e(x(tau-))` --
  and the output becomes one point evaluation at the end of `h` over those states, so the adjoint is one backward
  problem per output seeded with `dh/d(q, s)` and the forward sensitivities those of a point output. OCFESLV's
  capture of evolution-direction reductions in outputs (`f(cap)` at completion) is the counterpart on the
  collocation side. Under MoL a spatial reduction nested inside a nonlinearity (`(int u dz)^2`) becomes a scalar
  part composed at every `t` -- the collocation's root-level assembly, linear by construction, sums such an
  integrand per element (a third of the value on three elements, on OCFESLV's monolithic path too: an open item
  there, reproducer attached).
* **DAESLV: the collocation options declared once.** `DAESLV::Options` derives from `OCFEBASE::Options` as
  `OCFESLV::Options` does, so `options.INTERFACE`, `.AUDIT`, `.DETERMINACY`, `.FATAL` are the same path on both
  solvers; `options.MOL` keeps the method of lines' own group (`NEUMANN_NODE`, `REDUCE_HIGH_INDEX`). The second copy
  of the groups under `MOL` is gone (`options.MOL.INTERFACE.*` -> `options.INTERFACE.*`, `options.MOL.IMPOSITION.*`
  -> `options.MOL.*`). `OUTPUT` stays the integrator's group. In Python `DAESLV.Options` is a `CollocationOptions`,
  `DAESLV.MOLOptions` the MoL group, the `NeumannNode` enumerators on `DAESLV`.
* **OCFESLV: sensitivities by name.** `solve_fsens` / `solve_asens` take the inputs and constants as `solve` does,
  every given input a control of the solve, and leave the gradient in `val_function_gradient()`: one row per
  declared DOF, one column per output -- the layout DAESLV and ODESLV give, so the three solvers read the same.
  Bound in Python. OCFESLV had only the registry surface, `sens_jacobian()` over the control vector.
* **Python: a model distributed in space on DAESLV.** `DAESLV.semidiscrete_states()` and
  `DAESLV.semidiscrete_nodes(field)` expose the semi-discrete system's unknowns and a field's nodes, mesh order,
  as `(column, {direction: coordinate})`; everything else a PDE model needs was bound since the by-name solve. New
  `check_daeslv_mol.py`: a rod on DAESLV and on OCFESLV from one declaration, a time-varying inlet with one
  gradient row per declared DOF, a transition against closed forms, the semi-discrete layout.
* **OCFEMOL: the residual's Jacobian coloured on the equation rows alone.** The coloured pattern held the output
  rows too, for OCFESLV's single sweep; under MoL an output distributed over space is one row touching every
  node of its fields, so every such column needed its own colour: 121 for a residual whose rows need far fewer,
  paid at every Jacobian and every sensitivity or backward residual. OCFEMOL sweeps the residual and the outputs
  separately and now colours them separately. PSA column: `solve_fsens` 2.6 s -> 1.95 s, `solve_asens` 3.9 s -> 2.7 s.
* **DAESLV: `ic_sens` and `ic_adjoint` factorise as the IC solve does.** One `_ic_factor` builds the IC system's
  Jacobian in the pattern the IC solve uses, transposed for the adjoint, on KLU under `LINSOL = SPARSE` and dense
  otherwise; the sensitivity seeds at t0 and the sensitivities across a stop, forward and adjoint, were dense
  whatever `LINSOL` said.
* **DAESLV: the consistent initialisation follows the integrator's linear algebra.** With `LINSOL = SPARSE` the IC
  solve (and its polish) factorises with KLU; it factorised dense whatever `LINSOL` said. `DAESLV_mol`, whose leg 1
  drives a 12 000-unknown space-time system through DAESLV: 514 s -> 48 s in the sandbox.
* **DAESLV: IDACalcIC retired.** The KINSOL consistent initialisation is the one corrector, at t0 and at every
  stage boundary; its `(y0, y'0)` now starts IDAS (the derivative was computed and discarded, with IDACalcIC
  repairing it after the solve by default). Gone with it: `INITIAL.SOLVER`, `CALCIC`, `TOUT1`, `NMAXIC`,
  `MAXJACIC`, `MAXNLIIC`, `LSOFFIC`, `STEPTOLIC` and DAESLV's unused `BYPASS_EXPLICIT` (ODESLV keeps its own).
  `stats.numIC` counts the t0 solve. Results move at tolerance level: IDAS starts from KINSOL's point.
* **Quadratures accumulate over the horizon** in ODESLV and DAESLV, with no per-stage reset; an integral over t is
  the quadrature's value at the end. Fixes DAESLV failing at a stage boundary without a state restart (a large
  integrand against ATOL) and integral outputs that could come out wrong on a multi-stage horizon.
* **MoL on DAESLV (6a, 6b):** a model distributed in space declared through FFModel, solved by DAESLV with inputs,
  outputs and sensitivities by name; `tools/ci/build_deps.sh` patches SUNDIALS IDAS's order decision (an exactly
  constant state with an error-controlled quadrature stalled it at order 1).

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

* Requires MC++ 5.1 (`pymcpp~=5.1.0`): its `FFSubgraph` slot guard, without which the embedded solvers segfault
  when the model's own variables are the operation's inputs, and its cone-restricted `SFAD`/`SBAD`, which made
  CRONOS's own row-by-row differentiation (`_sfad_by_row`) redundant -- retired, the seven call sites going to
  `SFAD` directly, with the setup times equal within noise.

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
