# Changes

## Unreleased

### DAESLV: models distributed in space by the method of lines

* **DAESLV takes a model distributed in space**, declared through FFModel exactly as for OCFESLV, and integrates it by
  the method of lines: the declared model is validated as by any solver, semi-discretised internally (OCFEMOL), and
  the semi-discrete DAE becomes the working model -- one solver, the declared model staying the user's. What MoL does
  not take yet is refused by name at setup.
* **OCFEMOL, the semi-discretisation (MoL phases 1-5).** It shares the collocation base with OCFESLV (OCFEBASE, below):
  the layout, the interface plan, the DOF audit and the sparse Jacobian. The model is read along `t` as
  `M(x) d_t x = R` (`FFModel::peel_evolution`, the per-row peel hoisted from BASE_DE); the mass matrix comes from the row
  provenance with the node-level DOF audit and the INITIAL rows as a stream of scalar expressions; the structural index
  along `t` is checked on the full model, and an index-2 PDE is reduced to index 1 automatically under
  `REDUCE_HIGH_INDEX`; a model with no spatial direction (an ODE, a DAE) is the DAE itself. The residual and its sparse
  Jacobian (compressed forward AD over the base's colouring) enter the integrator's DAG as external operations
  (`FFMOLRES` / `FFGradMOLRES`); algebraic components are kept out of IDAS's error test; a field without a continuity
  claim follows the fold rule. Gates: DAESLV_molsetup, DAESLV_molgate (the heated rod and PDE5's 2-D heat equation
  against OCFESLV and exact solutions), DAESLV_molbench, DAESLV_molsens (the PSA column).
* **DAESLV: the collocation options declared once.** `DAESLV::Options` derives from `OCFEBASE::Options` as
  `OCFESLV::Options` does, so `options.INTERFACE`, `.AUDIT`, `.DETERMINACY`, `.FATAL` are the same path on both
  solvers; `options.MOL` keeps the method of lines' own group (`NEUMANN_NODE`, `REDUCE_HIGH_INDEX`). The second copy
  of the groups under `MOL` is gone (`options.MOL.INTERFACE.*` -> `options.INTERFACE.*`, `options.MOL.IMPOSITION.*`
  -> `options.MOL.*`). `OUTPUT` stays the integrator's group. In Python `DAESLV.Options` is a `CollocationOptions`,
  `DAESLV.MOLOptions` the MoL group, the `NeumannNode` enumerators on `DAESLV`.
* **MoL on DAESLV (6a, 6b):** a model distributed in space declared through FFModel, solved by DAESLV with inputs,
  outputs and sensitivities by name; `tools/ci/build_deps.sh` patches SUNDIALS IDAS's order decision (an exactly
  constant state with an error-controlled quadrature stalled it at order 1).
* **MoL on DAESLV (6c): time-varying inputs.** Each spatial slot of a declared input over `t` becomes an input over
  `t` of the semi-discrete DAE with the declared collocation; declared DOFs and working parameters are matched as
  a bijection, so by-name values and gradient rows read as the declared DOFs. `ic_sens` now differentiates through
  the placeholders binding a time-varying input at `t0` (an INITIAL row `x(0) = u(t0)` had no sensitivity to `u`,
  in ODESLV and DAESLV alike), and such a row no longer counts as explicit for the bypass.
* **MoL on DAESLV (6c): transitions.** A declared transition is applied at its state's differential nodes, each
  component composed pointwise there (`OCFEMOL::transition_expressions`, the INITIAL rows' form: the integrator
  peels these rows, so no vectorised operation); the integrator's stop machinery re-solves the algebraic nodes and
  carries the sensitivities across. `FFMOLINI`/`FFGradMOLINI`, the vectorised Ini stream, retired for that reason.
* **MoL: trajectory outputs.** An output distributed over `t` is, for each spatial row of its spatial part, an output
  distributed over `t` of the semi-discrete model: one value per stage time its mask selects, as the integrators
  take such an output; distributed over space too, one trajectory per spatial node, space-major.
* **MoL: transitions reading an input.** An input uniform in space -- a parameter, or an input over `t` only -- is
  bound to its semi-discrete slot when the transition is composed at the field's nodes; an input distributed in
  space is refused by name.
* **MoL: the consistent initialisation starts from the declared profile at every node.** An INITIAL row at an
  algebraic node (a Dirichlet node, or by default a Neumann node) is set aside, the node being held by its
  boundary row; `OCFEMOL::init_seed` skipped those rows, so the node started from its reference value, and a
  boundary row with several roots could be met on the wrong one. On the heated rod with `D(u) = D0 (1 + beta u)`
  the flux row was met at `1 + beta u < 0`, a negative diffusivity: `u(2,1) = -1.02` against OCFESLV's `2.79`.
  The set-aside rows now seed their nodes too; DAESLV reports, at `DISPLAY_LEVEL >= 1`, how far the consistent
  initialisation moved such a node from its declared value (`initial_consistency_defect`, unread before).
* **MoL: the structural index of the semi-discrete system is checked at setup** (the MoL design's step 6, missing
  until now): a maximum matching of the algebraic rows to the algebraic unknowns on `Rc`'s dependency pattern. Full:
  index 1. Deficient: index >= 2, refused by name -- an interface multiplier minted by the plan with no algebraic row,
  a constraint row reading differential unknowns only. Under `IC_TRACE` every value-continuity multiplier has that
  form (the Hessenberg index-2 pattern; KINSOL's IC Jacobian singular); under `IC_STRONG` the replaced copy is
  algebraic and pinned, index 1 -- verified on every MoL driver of the corpus. `semi_discrete().n_alg_rank` reports it.
* **MoL: null-space elimination of the plan's multipliers that raise the index** (phase 9, step 2). A multiplier the
  structural check finds with no algebraic row -- entering exactly two differential rows with opposite equal weights,
  its constraint a linear row `x_L - x_R = 0` on the two copies of a node -- is eliminated: the copies merged, their
  rows summed (the multiplier cancels), the constraint row and the multiplier column dropped. Exact and drift-free;
  applied as a layer, the base keeping its full layout and the integrator seeing the reduced unknowns. Keyed on the
  finding, not the imposition: every value-continuity multiplier under `IC_TRACE`, and any `IC_STRONG` claim that
  would land the same way; kept-explicit multipliers pinned by algebraic rows are left alone. On the rod under
  `IC_TRACE` the merged operator's spectrum is `IC_STRONG`'s to every digit, and so are the solution and its
  sensitivities. Anything else -- a nonlinear constraint, a multiplier joining more than two copies -- is refused by
  name. Gate `DAESLV_molindex`; `DAESLV_molsetup` G2. The copies a multiplier joins are named by its CONSTRAINT, not by
  its rows' masses: a balance row may carry several mass entries (the PSA column's c1 balance carries `F d_t q1`), and
  only the constrained unknowns are merged. On that column under `IC_TRACE` (six value multipliers, five fields, a
  nonlinear isotherm) the result equals `IC_STRONG` to 4e-11 on values and 4e-10 on forward and adjoint gradients.

### DAESLV: initialisation, sensitivities and adjoints

* **DAESLV: IDACalcIC retired.** The KINSOL consistent initialisation is the one corrector, at t0 and at every
  stage boundary; its `(y0, y'0)` now starts IDAS (the derivative was computed and discarded, with IDACalcIC
  repairing it after the solve by default). Gone with it: `INITIAL.SOLVER`, `CALCIC`, `TOUT1`, `NMAXIC`,
  `MAXJACIC`, `MAXNLIIC`, `LSOFFIC`, `STEPTOLIC` and DAESLV's unused `BYPASS_EXPLICIT` (ODESLV keeps its own).
  `stats.numIC` counts the t0 solve. Results move at tolerance level: IDAS starts from KINSOL's point.
* **DAESLV: the consistent initialisation converges to the integrator's tolerance.** KINSOL's point met the
  residual tolerance but pinned an ill-conditioned algebraic unknown only to ~1e-12, which IDAS error-tests at
  `ATOL`; a polish pass with KINSOL's scaled-step criterion in the error-weighted norm `1/(RTOL|y|+ATOL)` now ends
  every IC solve, at t0 and at each stop (`BASE_KINSOL::Knobs::RTOL/ATOL/WSTEPTOL`). This is what IDACalcIC did.
* **DAESLV: the consistent initialisation follows the integrator's linear algebra.** With `LINSOL = SPARSE` the IC
  solve (and its polish) factorises with KLU; it factorised dense whatever `LINSOL` said. `DAESLV_mol`, whose leg 1
  drives a 12 000-unknown space-time system through DAESLV: 514 s -> 48 s in the sandbox.
* **DAESLV: `ic_sens` and `ic_adjoint` factorise as the IC solve does.** One `_ic_factor` builds the IC system's
  Jacobian in the pattern the IC solve uses, transposed for the adjoint, on KLU under `LINSOL = SPARSE` and dense
  otherwise; the sensitivity seeds at t0 and the sensitivities across a stop, forward and adjoint, were dense
  whatever `LINSOL` said.
* **BASE_KINSOL: an automatic Newton-step bound** -- unbounded under the line search, KINSOL's own without it.
* **DAESLV adjoint: one backward problem for every adjoined output.** The outputs' multipliers are stacked into
  one IDA instance with one step sequence and one block-diagonal matrix (KLU whenever there is more than one
  block or the forward factorises sparse), so the forward rows and Jacobian are evaluated once per backward
  residual rather than once per output. The backward Jacobian is analytic (it was built by difference quotients),
  and the quadrature sensitivity integrands read the memoised evaluation. The PSA column at 480 states, one
  parameter, two outputs: `solve_fsens` 4.7 s -> 2.6 s, `solve_asens` 6.1 s -> 3.9 s, against 0.6 s for the solve.
* **DAESLV adjoint: a transition reading a time-varying input.** `ic_adjoint` lacked the placeholder term
  `-(dF/dq . dq/dp)^T mu` that `ic_sens` has had since 6c, so the adjoint gradient with respect to such an input
  was zero (`x(tau+) = x(tau-) + g`: forward `dJ/dg = (1, 0)`, adjoint `(0, 0)`; ODESLV was right). Pointwise and
  under MoL alike.
* **ODESLV adjoint: an INITIAL row reading a time-varying input** (`x(0) = u(t0)`) had a zero adjoint gradient in `u`:
  ODESLV closed its adjoint at t0 by differentiating `f0 + lambda . IC(p)` symbolically, and `u(t0)` reaches the IC
  expressions through a placeholder that derivative does not see. When the start system carries input placeholders,
  the initial-condition chain is now taken as DAESLV takes it -- `(lambda + df0/dx)(t0) . dx0/dp`, `dx0/dp` from
  `ic_sens` on the start system -- the symbolic route unchanged otherwise. DAESLV_moltime (C) now asserts ODESLV's
  adjoint (it reported it).
* **DAESLV adjoint: one backward matrix and linear solver per backward problem, all freed** -- all but the last leaked
  with two or more adjoined outputs (valgrind on DAESLV_quadacc). Backward problems honour the step controls.

### Outputs

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
* **Quadratures accumulate over the horizon** in ODESLV and DAESLV, with no per-stage reset; an integral over t is
  the quadrature's value at the end. Fixes DAESLV failing at a stage boundary without a state restart (a large
  integrand against ATOL) and integral outputs that could come out wrong on a multi-stage horizon.
* **OCFESLV: a spatial reduction nested inside a nonlinearity of an output is lifted into a hidden algebraic state.**
  The collocation evaluates an output element by element and sums at the root: exact for anything linear in a
  reduction, wrong for `(int u dz)^2` -- a third of the value on three equal elements, for every ratio of spatial
  integrals a user writes (a bed loading fraction). Mirroring OCFEMOL's spatial parts under the method of lines and
  BASE_DE's lift of nonlinear outputs, `OCFESLV::_lift_nested_reductions` (behind the new `FFModel::_on_declaration`
  hook, before any check processes the declared rows) rewrites the DECLARATION: each distinct spatial reduction R an
  output is not affine in (every dE/dR free of reductions) becomes a state J over the directions R does not consume,
  pinned by the algebraic row `J - R = 0` at every node (linear: exact), with a reference of the fields' times the
  integrated measure (a zero reference made the scaled Newton system singular), and the output reads J. Reductions
  consuming the evolution direction are left to FFModel's capture; a reduction consuming every direction of its
  fields is refused by name (a scalar hidden state, not yet). Gate OCFE_PDE41: the rod's `int J^2 dt`, `J(T)^2` and
  a ratio of spatial integrals to 1e-8 against closed forms; the PSA column's loading fraction, OCFESLV = DAESLV
  under MoL to 5e-8.
* **A row integrating over a direction it is not distributed over counted the field-free part of its integrand on
  ONE element** (`J - int (1 + u) dz = 0` over t pinned J to `int u dz + 1/n_el`). OCFEBASE's section 4a assembles
  such a row element by element and subtracts `(n_el - 1) g` for its non-integrated part g, which it evaluated with
  the integrated FIELDS zeroed -- leaving `int 1 dz` (or `int p dz`) standing on the one element used. g is now the
  residual with its FFIntegral NODES over the free directions replaced by 0 (one substitution at setup, exact for a
  residual linear in its integrals; a constant g broadcast). It reached outputs through FFModel's capture, which
  materialises a point-evaluated `int (1 + u) dz` as exactly such a row; integrated over t too it was exact (no g
  correction there). OCFE_PDE41 (C): `int (1 + u) dz`, `int (u + p) dz` and a ratio with a constant inside, at T, to
  1e-8 on three elements. The corpus drivers with integral rows (OCFE_scalar_*, OCFE_capture_nested, OCFE_byname,
  OCFE_hybrid_adequacy) are unchanged: none has a field-free part under such an integral.

### OCFESLV and the collocation base

* **OCFEBASE, the collocation base, in its own header `src/ocfebase.hpp`** -- the OCVar environment and evaluation
  kernels, the interface plan (frozen realisations, redundancy detection, trace projection, weak SAT terms, plan
  packing), the determinacy verdict, the collocated layout, the sparse Jacobian and the output evaluation
  (`_eval_fct`), moved out of OCFESLV and shared with OCFEMOL; the evolution window (ALL / ONE / NONE) replaces the
  marching hooks. Four dead OCFESLV internals removed.
* **Options split:** `OCFEBASE::Options` holds `INTERFACE`, `OUTPUT`, `AUDIT`, `DETERMINACY` and `FATAL`;
  `OCFESLV::Options` derives from it and adds `SOLVE`, `SENSITIVITY` and `HOMOTOPY`. The spelling at the user's side
  (`options.INTERFACE.IMPOSITION`, ...) is unchanged.
* **OCFESLV: `solve` by name**, as ODESLV and DAESLV (`FFModel::resolve_values` shared; coordinate generators;
  `input_dofs` / `input_ndof`).
* **OCFESLV: sensitivities by name.** `solve_fsens` / `solve_asens` take the inputs and constants as `solve` does,
  every given input a control of the solve, and leave the gradient in `val_function_gradient()`: one row per
  declared DOF, one column per output -- the layout DAESLV and ODESLV give, so the three solvers read the same.
  Bound in Python. OCFESLV had only the registry surface, `sens_jacobian()` over the control vector.
* **An energy-stable SAT form for IC_WEAK** (`INTERFACE.SAT_FORM = SAT_ENERGY`), alongside the legacy one, which stays
  OCFESLV's default. The legacy SAT gives the semi-discrete operator a growing mode quadratic in `SAT_SIGMA0` -- the
  value-jump penalty enters the auxiliary flux and the flux-jump penalty multiplies it (`+9.2e4` at 10 on the rod) --
  damped on OCFESLV by Radau's L-stability, fatal to an integrator. The energy form is the summation-by-parts one:
  the value jump lifted into the auxiliary's LINK row and the flux jump into the PDE row (`SAT_FLUX`: alternating by
  default, central available), an upwind term for a constant-coefficient convective part, an optional interior
  penalty (`SAT_PENALTY`, default 0). Every strength, legacy or energy, is computed in ONE function,
  `OCFEBASE::_sat_strength`; the frozen interface plan records the form, flux and penalty. Under the method of lines
  IC_WEAK takes the energy form whatever `SAT_FORM` says. Measured: stable for any penalty and either boundary
  arrangement; on pure diffusion one order better than IC_STRONG at even degree; on convection-diffusion behind it
  (2.6x at Peclet 10 on the rod, 20-25x on the PSA column) -- the open question. Gate `DAESLV_molsat`; G2 updated.
* **SAT_ENERGY: the upwind term of a convective part was generated only when SAT_PENALTY != 0** -- it sat in the
  interior-penalty block, switched off with the penalty when its default became 0 -- so the delivered energy form had
  NO upwind term: a growing mode at Peclet 100 (+0.46 in the semi-discrete spectrum; u(20,1) = -1240 on the rod under
  MoL). The block now runs for every SAT_ENERGY model and gates only the penalty. Peclet 10 on 6 elements: 4.75e-5 ->
  3.78e-5 (IC_STRONG 1.51e-5). Gate DAESLV_molsat (B2): Peclet 100 to t = 20 settles to u = 1. Found by the
  convective-gap probe, which also located the remaining gap AT THE INTERFACES (one element: every variant equal to the
  digit) and its mechanism: after order reduction the PDE row reads U q, q the auxiliary whose LINK row receives the
  value-jump lifting -- so the diffusive flux choice also sets a convective interface term.
* **SAT_ENERGY: the central flux is refused on a field with a convective part** (it diverges under MoL: 5e28 at 12
  elements, Peclet 10, with or without the upwind term), named at setup; taken on pure diffusion. Gate DAESLV_molsat
  (D). CLOSES the convective-gap probe (9-10 Nov 2026), measured on the rod and the PSA column against fine IC_STRONG
  references, degree 4: the gap is entirely at the INTERFACES (one element: every variant equal to the digit); it is
  not specific to the column (without adsorption the column's c1 is plain convection-diffusion, IC_WEAK's error stays
  ~4e-6 at 12 elements while IC_STRONG's falls to 4e-9 -- the energy form's interface error is roughly constant across
  problems, IC_STRONG exploits each problem's smoothness); the diffusive side is not the lever (SAT_PENALTY 0..100
  makes it slightly worse); the convective side is, but only so far -- the upwind term at 2x and 4x the textbook
  |U|/w narrows the column's gap from 19x to 7x, at 8x and 16x the rate collapses (0.6) and the error grows; cancelling
  the lifting's convective effect in the PDE row (the row reads U q) makes it slightly worse. The textbook |U|/w is
  kept. Conclusion: the limit is structural to the weak interface treatment of convection here; for convective models
  under MoL, IC_STRONG and IC_TRACE remain the recommended impositions.
* **OCPlan:** the aux-implied rule only where a face has one tangential direction -- fixes OCFESLV on `(t,x,y)`
  models when a time element has more nodes than a spatial one; the tangential-derivative implication read from the
  plan's own decisions; the definition-row lookup supplies every leaf of a LINK row.

### Index reduction

* **High-index reduction: a declared initial condition IMPLIED by the materialised hidden constraints is set aside.**
  With `REDUCE.HIDDEN_IC` the reduction materialises its hidden constraint levels as INITIAL rows at the initial
  point; a declared row they already imply over-determined it (OCFE_PDE20b, `x' = y, 0 = x - a(t)` with `x(0) = a(0)`
  declared: 13 rows for 12 unknowns, refused by the DOF audit -- the gate printed BASELINE since Stage 2). FFModel now
  matches the INITIAL rows to the differential states, the materialised levels first, and sets aside the declared rows
  left unmatched (`initial_set_aside()`, `initial_data().set_aside`) -- for every solver. Their consistency with the
  solution is checked: `OCFESLV::initial_consistency_defect()`, and both OCFESLV and DAESLV warn at DISPLAY_LEVEL >= 1
  when a declared value contradicts the hidden constraints (which hold). OCFE_PDE20b now passes (REDUCED:PASS, exact
  to 7e-16); new gate OCFE_PDE20h -- consistent, inconsistent (defect 0.5 reported), and alongside free data (kept).
  COVERAGE: a declared row covering MORE nodes than the level implying it (OCFE_PDE39b: IC_c on every z node, the level
  from an EOS row absent at z = LB) is NARROWED to the nodes the level leaves free instead of set aside whole (which
  removed one node too many: 215 rows for 216 unknowns); any other mismatch is left to the DOF audit. OCFE_PDE20e (B)
  now gates this convention (it gated the opposite: the declared IC rejected as a surplus).

### Python and notebooks

* **Python: one by-name convention for every solver** -- a value is a number, a list, or a function of the
  coordinates; `OCFESLV.solve` by name; `input_coordinates`.
* **Python: a model distributed in space on DAESLV.** `DAESLV.semidiscrete_states()` and
  `DAESLV.semidiscrete_nodes(field)` expose the semi-discrete system's unknowns and a field's nodes, mesh order,
  as `(column, {direction: coordinate})`; everything else a PDE model needs was bound since the by-name solve. New
  `check_daeslv_mol.py`: a rod on DAESLV and on OCFESLV from one declaration, a time-varying inlet with one
  gradient row per declared DOF, a transition against closed forms, the semi-discrete layout.
* **Notebooks: `MOL_tutorial` and `PSA_tutorial`.** The OCFESLV tutorial's heated rod on DAESLV by the method of
  lines, one declaration for both solvers: outputs, profiles, OCFESLV's spectral convergence in its time order,
  gradients by name in the shared layout. A non-isothermal two-component adsorption column, by MoL or by space-time
  collocation behind one switch: breakthrough, thermal front, sensitivities to the affinity, a cross-check.

### Performance

* **OCFEMOL: the residual's Jacobian coloured on the equation rows alone.** The coloured pattern held the output
  rows too, for OCFESLV's single sweep; under MoL an output distributed over space is one row touching every
  node of its fields, so every such column needed its own colour: 121 for a residual whose rows need far fewer,
  paid at every Jacobian and every sensitivity or backward residual. OCFEMOL sweeps the residual and the outputs
  separately and now colours them separately. PSA column: `solve_fsens` 2.6 s -> 1.95 s, `solve_asens` 3.9 s -> 2.7 s.
* **OCFEBASE: `_defining_rows()` computed once per plan build**, not per node.

### Build, tests and tools

* **Build and docs follow `ocfebase.hpp`:** install lists, the options generator (`src/interface/gen_options.py`,
  which formats with MC++'s `.clang-format` when `MCPP_ROOT` is set), the Doxygen input and the environment table.
  `check_daeslv_mol.py` runs with the other Python checks under `make check`; `check_ffodeslv` / `check_ffdaeslv`: the
  coordinate case no longer shadows the solver.
* **`tools/ci/snapshot_corpus.sh`:** a driver without a baseline entry is NEW, recorded on first sight and gated from
  the next check; an advisory names the drivers whose output prints a FAIL verdict line; the concurrency guard names
  the processes it saw; DAESLV and ODESLV timing families; clocks suppressed under `CRONOS_SNAPSHOT`.
* **New gates:** DAESLV_molsetup, DAESLV_molpeel, DAESLV_molgate, DAESLV_molbench, DAESLV_mol, DAESLV_molpde,
  DAESLV_molpsa, DAESLV_molsens, DAESLV_moltime, DAESLV_molindex, DAESLV_molsat; OCFE_PDE20h, OCFE_PDE41;
  ODESLV_evalsum's NONLINEAR case.

### Known gaps

* **`SAT_ENERGY` on convective models:** stable, but behind `IC_STRONG` in accuracy (3x on a Peclet-10 rod, ~19x on
  the PSA column at 12 elements of degree 4) -- structural to the weak interface treatment of convection here; for
  convective models under MoL, `IC_STRONG` and `IC_TRACE` are recommended.
* **MoL: transitions reading an input distributed in space** are refused by name.
* **OCFESLV: a nested spatial reduction consuming every direction of its fields** (a scalar hidden state) is refused
  by name.

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
