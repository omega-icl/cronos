# %% [markdown]
# # CRONOS Python bindings -- `DAESLV` check
#
# Exercises `cronos.DAESLV` (CVODES): solves with and without forward and adjoint
# sensitivities, inputs given positionally and by name, transitions, distributed
# controls, results and statistics, options. Values are checked against the C++
# driver `ODESLV_test1` on the same model and against closed forms. Run as a script,
# it exits with status 1 on any failure.
#
# Requires `pymcpp` and `cronos` on the Python path (`CRONOS_PYPATH`).

# %%
import math
import os
import sys

for d in reversed(os.environ.get("CRONOS_PYPATH", "").split(":")):
    if d and d not in sys.path:
        sys.path.insert(0, d)

from pymcpp import FFGraph, FFPartial, FFIntegral, FFEval
from cronos import FFDom, FFModel, DAESLV

npass = nfail = 0


def check(cond, what):
    """Record and print one check."""
    global npass, nfail
    print("  %s  %s" % ("PASS" if cond else "FAIL", what))
    if cond:
        npass += 1
    else:
        nfail += 1


def close(a, b, rtol):
    """Element-wise relative closeness of two nested lists of floats."""
    if isinstance(a, (list, tuple)):
        return len(a) == len(b) and all(close(x, y, rtol) for x, y in zip(a, b))
    return abs(a - b) <= rtol * max(1., abs(b))


OpP, OpI, OpE = FFPartial(), FFIntegral(), FFEval()
INT = FFDom.ALL - FFDom.LB

# %% [markdown]
# ## A. The model of `ODESLV_test1`, against the C++ driver
#
# Reference values printed by `ODESLV_test1` (7 significant digits): outputs
# `f = [10.05033, 0.8200303]`; forward-sensitivity gradient
# `fp = [[-0.7158895, 1.974535], [-0.001896783, 0.0007985682]]`.

# %%
G = FFGraph()
t = G.add_var("t")
x0, x1 = G.add_var("x0(t)"), G.add_var("x1(t)")
p0, p1 = G.add_var("p0"), G.add_var("p1")
A = DAESLV(G)
A.options.SOLVE.VERBOSE = 0
A.add_domain(t, FFDom(0., 10., 4, FFDom.LGR, 4))
A.add_state(x0, [t])
A.add_state(x1, [t])
A.add_input(p0)
A.add_input(p1)
A.set_evolution_domain(t)
A.update_ref(x0, 1.2)
A.update_ref(x1, 1.1)
io = FFModel.EqnOptions(FFModel.EqnRole.INTERIOR)
ii = FFModel.EqnOptions(FFModel.EqnRole.INITIAL)
A.add_equation(OpP(x0, t) - p0 * x0 * (1. - x1), [t], [INT], io)
A.add_equation(OpP(x1, t) - p0 * x1 * (x0 - 1.), [t], [INT], io)
A.add_equation(x0 - 1.2, [t], [FFDom.LB], ii)
A.add_equation(x1 - (1.1 + 0.01 * p1), [t], [FFDom.LB], ii)
A.add_output(OpI(x1, t))
A.add_output(OpE(x0 * x1, t, 10.))
o = A.options
o.SOLVE.INTMETH = DAESLV.Options.MSBDF
o.SOLVE.NLINSOL = DAESLV.Options.NEWTON
o.SOLVE.LINSOL = DAESLV.Options.SPARSE
o.SOLVE.FSACORR = DAESLV.Options.STAGGERED
o.SOLVE.NMAX = 2000
o.SOLVE.ATOL = o.SOLVE.ATOLB = o.SOLVE.ATOLS = 1e-9
o.SOLVE.RTOL = o.SOLVE.RTOLB = o.SOLVE.RTOLS = 1e-9
o.SOLVE.QERR = o.SOLVE.QERRS = True
o.SOLVE.ASACHKPT = 2000
o.SOLVE.VERBOSE = 0
check(A.setup(), "setup()  (%s)" % (A.extract_error() or "ok"))
check(A.np() == 2 and A.nf() == 2 and A.nx() == 2, "sizes: np = nf = nx = 2")

F_REF = [10.05033, 0.8200303]
FP_REF = [[-0.7158895, 1.974535], [-0.001896783, 0.0007985682]]
p = [2.96, 3.]
check(A.solve(p) == DAESLV.NORMAL, "solve(p) -> NORMAL")
f = A.val_function()
print("  f =", f)
check(close(f, F_REF, 2e-6), "outputs == C++ driver (7 digits)")

check(A.solve_fsens(p) == DAESLV.NORMAL, "solve_fsens(p) -> NORMAL")
gf = A.val_function_gradient()
print("  forward gradient =", gf)
gf_t = [list(r) for r in zip(*gf)]
match = "as printed" if close(gf, FP_REF, 5e-6) else ("transposed" if close(gf_t, FP_REF, 5e-6) else None)
check(match == "as printed", "forward gradient == C++ fp, one row per direction (p0, p1), one column per output")

check(A.solve_asens(p) == DAESLV.NORMAL, "solve_asens(p) -> NORMAL")
ga = A.val_function_gradient()
print("  adjoint gradient =", ga)
check(close(ga, gf, 1e-6), "adjoint gradient == forward gradient (1e-6)")
st = A.stats_solve
print("  stats_solve: steps %d, RHS %d, JAC %d | stats_fsens: steps %d | stats_asens: steps %d"
      % (st.numSteps, st.numRES, st.numJAC, A.stats_fsens.numSteps, A.stats_asens.numSteps))
check(st.numSteps > 0 and A.stats_fsens.numSteps > 0 and A.stats_asens.numSteps > 0,
      "statistics recorded separately for solve, solve_fsens and solve_asens")

# %% [markdown]
# ## B. A closed form, inputs by name, recorded results
#
# x' = -p x, x(0) = 1: x(1) = exp(-p) and dx(1)/dp = -exp(-p).

# %%
G2 = FFGraph()
t, x, p = G2.add_var("t"), G2.add_var("x(t)"), G2.add_var("p")
B = DAESLV(G2)
B.options.SOLVE.VERBOSE = 0
B.add_domain(t, FFDom(0., 1., 2, FFDom.LGR, 3))
B.set_evolution_domain(t)
B.add_state(x, [t], ref=1.)
B.add_input(p)
B.add_equation(OpP(x, t) + p * x, [t], [INT])
B.add_equation(x - 1., [t], [FFDom.LB], ii)   # INITIAL: AUTO would resolve LB alone to BOUNDARY
B.add_output(OpE(x, t, 1.))
B.options.SOLVE.RTOL = B.options.SOLVE.ATOL = B.options.SOLVE.RTOLS = B.options.SOLVE.ATOLS = 1e-10
B.options.SOLVE.RTOLB = B.options.SOLVE.ATOLB = 1e-10
B.options.OUTPUT.RESRECORD = 10
check(B.setup(), "setup()")
check(B.solve({p: [0.7]}) == DAESLV.NORMAL, "solve by name {p: [0.7]}")
check(close(B.val_function(), [math.exp(-0.7)], 1e-8), "x(1) == exp(-0.7)")
check(B.solve({p: 0.7}) == DAESLV.NORMAL and close(B.val_function(), [math.exp(-0.7)], 1e-8),
      "solve by name {p: 0.7}: a bare number for a one-value input")
check(B.solve_fsens({p: [0.7]}) == DAESLV.NORMAL and close(B.val_function_gradient(), [[-math.exp(-0.7)]], 1e-7),
      "forward: dx(1)/dp == -exp(-0.7)")
check(B.solve_asens({p: [0.7]}) == DAESLV.NORMAL and close(B.val_function_gradient(), [[-math.exp(-0.7)]], 1e-7),
      "adjoint: dx(1)/dp == -exp(-0.7)")
import numpy as np
B.solve({p: [0.7]})
r = B.results_solve
print("  results_solve: t %s, x %s, q %s" % (r.t.shape, r.x.shape, r.q.shape))
check(isinstance(r.x, np.ndarray) and r.x.shape == (len(r.t), 1) and r.q.shape == (len(r.t), 0),
      "results_solve: NumPy arrays t (n,), x (n, nx), q (n, nq)")
check(abs(r.t[-1] - 1.) < 1e-12 and np.allclose(r.x[:, 0], np.exp(-0.7 * r.t), rtol=1e-7),
      "results_solve: x(t) == exp(-0.7 t) along the trajectory")
B.solve_fsens({p: [0.7]})
rf = B.results_fsens
check(rf.xp.shape == (1, len(rf.t), 1) and np.allclose(rf.xp[0, :, 0], -rf.t * np.exp(-0.7 * rf.t), rtol=1e-6, atol=1e-9),
      "results_fsens: xp (nsen, n, nx) == -t exp(-0.7 t)")
B.solve_asens({p: [0.7]})
ra = B.results_asens
check(ra.l.shape == (1, len(ra.t), 1) and ra.qp.shape == (1, len(ra.t), 1) and np.all(np.diff(ra.t) >= 0),
      "results_asens: l (nf, n, nx), qp (nf, n, nsen), time non-decreasing (stage boundaries recorded twice)")
check(abs(ra.qp[0, 0, 0] + math.exp(-0.7)) < 1e-7 and abs(ra.qp[0, -1, 0]) < 1e-12,
      "results_asens: qp(0) == gradient -exp(-0.7), qp(1) == 0")
check(len(B.results_fsens.t) == len(rf.t), "results_fsens kept after solve_asens (separate records)")
try:
    B.solve({})
    check(False, "an incomplete by-name dict is refused")
except ValueError as e:
    check(True, "an incomplete by-name dict is refused: %s" % str(e)[:60])

# %% [markdown]
# ## C. A transition
#
# x' = -x, x(0) = 1, jump x(0.5^+) = x(0.5^-) + 1: x(1) = exp(-1) + exp(-0.5).

# %%
G3 = FFGraph()
t, x = G3.add_var("t"), G3.add_var("x(t)")
C = DAESLV(G3)
C.options.SOLVE.VERBOSE = 0
C.add_domain(t, FFDom(0., 1., 4, FFDom.LGR, 3))
C.set_evolution_domain(t)
C.add_state(x, [t], ref=1.)
C.add_equation(OpP(x, t) + x, [t], [INT])
C.add_equation(x - 1., [t], [FFDom.LB], ii)
C.add_transition(x + 1., x, t, .5)
C.add_output(OpE(x, t, 1.))
C.options.SOLVE.RTOL = C.options.SOLVE.ATOL = 1e-10
check(C.setup(), "setup() with a transition  (%s)" % (C.extract_error() or "ok"))
check(C.solve([]) == DAESLV.NORMAL and close(C.val_function(), [math.exp(-1.) + math.exp(-0.5)], 1e-8),
      "x(1) == exp(-1) + exp(-0.5)")

# %% [markdown]
# ## D. A distributed control
#
# x' = u(t), x(0) = 0, u piecewise constant on 4 elements of [0, 1]:
# x(1) = sum of u_e / 4 and dx(1)/du_e = 1/4. Values by DOF index are a list
# built over input_dofs; a callable is always a function of the coordinates.

# %%
G4 = FFGraph()
t, x, u = G4.add_var("t"), G4.add_var("x(t)"), G4.add_var("u(t)")
D = DAESLV(G4)
D.options.SOLVE.VERBOSE = 0
D.add_domain(t, FFDom(0., 1., 4, FFDom.LGR, 3))
D.set_evolution_domain(t)
D.add_state(x, [t], ref=0.)
D.add_input(u, [t], type=FFDom.LG, n_node=1, ref=1.)
D.add_equation(OpP(x, t) - u, [t], [INT])
D.add_equation(x, [t], [FFDom.LB], ii)
D.add_output(OpE(x, t, 1.))
D.register_control(u)
D.options.SOLVE.RTOL = D.options.SOLVE.ATOL = D.options.SOLVE.RTOLS = D.options.SOLVE.ATOLS = 1e-10
check(D.setup(), "setup() with a distributed control  (%s)" % (D.extract_error() or "ok"))
check(D.n_control_dof() == 4, "4 control DOFs (one per element)")
vals = [1. + d.element[t] for d in D.input_dofs(u)]          # u_e = 1, 2, 3, 4: values by DOF index
check(D.solve_fsens({u: vals}) == DAESLV.NORMAL and close(D.val_function(), [2.5], 1e-8),
      "x(1) == (1+2+3+4)/4 = 2.5  (values by DOF index, a list over input_dofs)")
g = D.val_function_gradient(u)
print("  dx(1)/du =", g)
check(close(g, [[.25]] * 4, 1e-7), "dx(1)/du_e == 1/4 for every element (one row per DOF)")
P = D.set_parameter_values(u, [1., 2., 3., 4.], [0.] * D.np())
check(close(D.get_parameter_values(u, P), [1., 2., 3., 4.], 0.), "set/get_parameter_values round trip")
mids = [c[t] for c in D.input_coordinates(u)]
check(close(mids, [.125, .375, .625, .875], 1e-14), "input_coordinates: one LG node per element, the midpoints")
check(D.solve({u: lambda coord: 1. + 4. * coord[t]}) == DAESLV.NORMAL and close(D.val_function(), [3.], 1e-8),
      "a callable is a function of the coordinates: u = 1 + 4t at the midpoints, x(1) == 3")
try:
    D.solve({u: lambda dof: 1. + dof.element[t]})
    check(False, "a DofIndex-style callable is refused")
except ValueError as e:
    check("input_dofs(" in str(e), "a DofIndex-style callable is refused, pointing to input_dofs")

# %% [markdown]
# ## E. Options

# %%
o = DAESLV.Options()
print("  defaults: INTMETH", o.SOLVE.INTMETH, "| LINSOL", o.SOLVE.LINSOL, "| RTOL", o.SOLVE.RTOL, "| REDUCE.ORDER", o.REDUCE.ORDER)
check(isinstance(o, FFModel.Options), "DAESLV.Options extends FFModel.Options")
print("  doc of RTOL:", DAESLV.Options.t_Solve.RTOL.__doc__)   # the nested GROUP type, not the property

# %% [markdown]
# ## blk_fct and the output index

# %%
_b = [A.blk_fct(k) for k in range(2)]
check(_b == [(0, 1), (1, 1)], "blk_fct(k) = (k, 1): an DAESLV output is ONE value (%s)" % str(_b))
try:
    A.blk_fct(2); check(False, "blk_fct beyond the last output raises IndexError")
except IndexError:
    check(True, "blk_fct beyond the last output raises IndexError")
_G6 = FFGraph(); _t6, _x6, _k6 = _G6.add_var("t"), _G6.add_var("x"), _G6.add_var("k")
_S6 = DAESLV(_G6)
_S6.options.SOLVE.VERBOSE = 0
_i0 = _S6.add_output(OpE(_x6, _t6, 1.)); _i1 = _S6.add_output(OpI(_x6, _t6))
check((_i0, _i1) == (0, 1), "add_output returns the output's index on an DAESLV too (%s)" % str((_i0, _i1)))

# %% [markdown]
# ## An output distributed over the evolution direction: one value per stage time

# %%
_G8 = FFGraph(); _t8, _x8, _k8 = _G8.add_var("t"), _G8.add_var("x"), _G8.add_var("k")
_S8 = DAESLV(_G8)
_S8.add_domain(_t8, FFDom(0., 1., 4, FFDom.LGR, 3)); _S8.set_evolution_domain(_t8)
_S8.add_state(_x8, [_t8], ref=1.); _S8.add_input(_k8, ref=0.8)
_S8.add_equation(OpP(_x8, _t8) + _k8 * _x8, [_t8], [FFDom.ALL - FFDom.LB], DAESLV.EqnOptions(DAESLV.EqnRole.INTERIOR))
_S8.add_equation(_x8 - 2. * _k8, [_t8], [FFDom.LB], DAESLV.EqnOptions(DAESLV.EqnRole.INITIAL))   # x(t) = 2 k exp(-k t)
_j0 = _S8.add_output(_x8, [_t8], masks=[FFDom.ALL])
_j1 = _S8.add_output(OpE(_x8, _t8, 1.))
_S8.options.SOLVE.VERBOSE = 0
check(_S8.setup(), "a model with an output distributed over t sets up (%s)" % (_S8.extract_error() or "ok"))
check((_j0, _j1) == (0, 1) and _S8.blk_fct(0) == (0, 5) and _S8.blk_fct(1) == (5, 1),
      "blk_fct: the distributed output takes one value per stage time (5), the point output the next (%s, %s)" % (_S8.blk_fct(0), _S8.blk_fct(1)))
_TS = np.array([0., .25, .5, .75, 1.])
_ex = 2 * .8 * np.exp(-.8 * _TS); _dex = 2 * np.exp(-.8 * _TS) * (1 - .8 * _TS)
_S8.solve_fsens([0.8]); _f = np.array(_S8.val_function()); _g = np.array(_S8.val_function_gradient()).ravel()
check(len(_f) == 6 and np.abs(_f[:5] - _ex).max() < 1e-6 and abs(_f[5] - _ex[-1]) < 1e-6,
      "the values at the 5 stage times = 2 k exp(-k t), and the point output after them")
check(len(_g) == 6 and np.abs(_g[:5] - _dex).max() < 1e-5, "forward gradient at the stage times = the closed form")
_S8.solve_asens([0.8]); _ga = np.array(_S8.val_function_gradient()).ravel()
check(len(_ga) == 6 and np.abs(_ga[:5] - _dex).max() < 1e-5,
      "adjoint gradient at the stage times = the closed form -- including the initial time (d x0 / dk = 2)")

# %% [markdown]
# ## G. A model distributed in space: the method of lines
#
# The heated rod of DAESLV_molpde: u(t, z) on z in [0, 1], u = 0 at z = 0, a unit flux at z = 1.
# DAESLV semi-discretises it at setup -- its working model is the semi-discrete DAE, one state per
# nodal unknown -- options.MOL holds the method of lines' own group, the collocation groups sit on options.

# %%
import cronos
G5 = FFGraph()
t, z, u = G5.add_var("t"), G5.add_var("z"), G5.add_var("u(t,z)")
R = DAESLV(G5)
R.options.SOLVE.VERBOSE = 0
R.options.DISPLAY_LEVEL = 0
R.options.SOLVE.RTOL, R.options.SOLVE.ATOL = 1e-10, 1e-12
R.options.OUTPUT.RESRECORD = 10
R.add_domain(t, FFDom(0., 1., 2, FFDom.LGR, 3))                    # the stages
R.add_domain(z, FFDom(0., 1., 3, FFDom.LGL, 5))
R.set_evolution_domain(t)
R.add_state(u, [t, z])
OpP = FFPartial()
TI, ZI = FFDom.ALL - FFDom.LB, FFDom.ALL - FFDom.LB - FFDom.UB
R.add_equation(OpP(u, t) - OpP(0.1 * OpP(u, z), z), [t, z], [TI, ZI], FFModel.EqnOptions(FFModel.EqnRole.INTERIOR))
R.add_equation(u, [t, z], [TI, FFDom.LB], FFModel.EqnOptions(FFModel.EqnRole.BOUNDARY))
R.add_equation(0.1 * OpP(u, z) - 1., [t, z], [TI, FFDom.UB], FFModel.EqnOptions(FFModel.EqnRole.BOUNDARY))
R.add_equation(u, [t, z], [FFDom.LB, FFDom.ALL], FFModel.EqnOptions(FFModel.EqnRole.INITIAL))

check(isinstance(R.options.MOL, DAESLV.MOLOptions) and isinstance(R.options, cronos.CollocationOptions),
      "options.MOL is a DAESLV.MOLOptions (the method of lines' own group); options itself a CollocationOptions")
check(R.options.INTERFACE.IMPOSITION == cronos.OCFESLV.Options().INTERFACE.IMPOSITION and hasattr(R.options, "AUDIT"),
      "the collocation groups (INTERFACE, AUDIT, ...) sit on options, at OCFESLV's path")
check(isinstance(cronos.OCFESLV.Options(), cronos.CollocationOptions) and isinstance(cronos.OCFESLV.Options(), FFModel.Options),
      "OCFESLV.Options derives from CollocationOptions (and so still from FFModel.Options)")
check(R.options.MOL.NEUMANN_NODE == DAESLV.NEUMANN_ALGEBRAIC,
      "options.MOL.NEUMANN_NODE defaults to NEUMANN_ALGEBRAIC")
check(R.setup() and R.solve([]) == DAESLV.NORMAL, "a model distributed in space sets up and solves on DAESLV")
xa = R.results_solve.x.copy()
check(xa.ndim == 2 and xa.shape[1] == 30, "results_solve: the trajectory of the 30 nodal unknowns (3 elements x 5 LGL in z, u and Dz_u)")
R.options.MOL.NEUMANN_NODE = DAESLV.NEUMANN_DIFFERENTIAL
check(R.setup() and R.solve([]) == DAESLV.NORMAL and not np.array_equal(R.results_solve.x[-1], xa[-1]),
      "options.MOL reaches the semi-discretisation (NEUMANN_DIFFERENTIAL changes the result)")
check(R.distributed_in_space(), "the declared model is still the PDE after setup")

# %% [markdown]
# ## Summary

# %%
print("\ncheck_daeslv: %d passed, %d failed -- %s" % (npass, nfail, "FAILURES" if nfail else "ALL PASS"))
if __name__ == "__main__" and "ipykernel" not in sys.modules and nfail:
    sys.exit(1)
