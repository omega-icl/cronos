# %% [markdown]
# # CRONOS Python bindings -- `ODESLV` check
#
# Exercises `cronos.ODESLV` (CVODES): solves with and without forward and adjoint
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
from cronos import FFDom, FFModel, ODESLV

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
A = ODESLV(G)
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
o.INTMETH = ODESLV.Options.MSBDF
o.NLINSOL = ODESLV.Options.NEWTON
o.LINSOL = ODESLV.Options.SPARSE
o.FSACORR = ODESLV.Options.STAGGERED
o.NMAX = 2000
o.ATOL = o.ATOLB = o.ATOLS = 1e-9
o.RTOL = o.RTOLB = o.RTOLS = 1e-9
o.QERR = o.QERRS = True
o.ASACHKPT = 2000
o.DISPLAY = 0
check(A.setup(), "setup()  (%s)" % (A.extract_error() or "ok"))
check(A.np() == 2 and A.nf() == 2 and A.nx() == 2, "sizes: np = nf = nx = 2")

F_REF = [10.05033, 0.8200303]
FP_REF = [[-0.7158895, 1.974535], [-0.001896783, 0.0007985682]]
p = [2.96, 3.]
check(A.solve(p) == ODESLV.NORMAL, "solve(p) -> NORMAL")
f = A.val_function()
print("  f =", f)
check(close(f, F_REF, 2e-6), "outputs == C++ driver (7 digits)")

check(A.solve_fsens(p) == ODESLV.NORMAL, "solve_fsens(p) -> NORMAL")
gf = A.val_function_gradient()
print("  forward gradient =", gf)
gf_t = [list(r) for r in zip(*gf)]
match = "as printed" if close(gf, FP_REF, 5e-6) else ("transposed" if close(gf_t, FP_REF, 5e-6) else None)
check(match == "as printed", "forward gradient == C++ fp, one row per direction (p0, p1), one column per output")

check(A.solve_asens(p) == ODESLV.NORMAL, "solve_asens(p) -> NORMAL")
ga = A.val_function_gradient()
print("  adjoint gradient =", ga)
check(close(ga, gf, 1e-6), "adjoint gradient == forward gradient (1e-6)")
st = A.stats_solve
print("  stats_solve: steps %d, RHS %d, JAC %d | stats_fsens: steps %d | stats_asens: steps %d"
      % (st.numSteps, st.numRHS, st.numJAC, A.stats_fsens.numSteps, A.stats_asens.numSteps))
check(st.numSteps > 0 and A.stats_fsens.numSteps > 0 and A.stats_asens.numSteps > 0,
      "statistics recorded separately for solve, solve_fsens and solve_asens")

# %% [markdown]
# ## B. A closed form, inputs by name, recorded results
#
# x' = -p x, x(0) = 1: x(1) = exp(-p) and dx(1)/dp = -exp(-p).

# %%
G2 = FFGraph()
t, x, p = G2.add_var("t"), G2.add_var("x(t)"), G2.add_var("p")
B = ODESLV(G2)
B.add_domain(t, FFDom(0., 1., 2, FFDom.LGR, 3))
B.set_evolution_domain(t)
B.add_state(x, [t], ref=1.)
B.add_input(p)
B.add_equation(OpP(x, t) + p * x, [t], [INT])
B.add_equation(x - 1., [t], [FFDom.LB], ii)   # INITIAL: AUTO would resolve LB alone to BOUNDARY
B.add_output(OpE(x, t, 1.))
B.options.RTOL = B.options.ATOL = B.options.RTOLS = B.options.ATOLS = 1e-10
B.options.RTOLB = B.options.ATOLB = 1e-10
B.options.RESRECORD = 10
check(B.setup(), "setup()")
check(B.solve({p: [0.7]}) == ODESLV.NORMAL, "solve by name {p: [0.7]}")
check(close(B.val_function(), [math.exp(-0.7)], 1e-8), "x(1) == exp(-0.7)")
check(B.solve_fsens({p: [0.7]}) == ODESLV.NORMAL and close(B.val_function_gradient(), [[-math.exp(-0.7)]], 1e-7),
      "forward: dx(1)/dp == -exp(-0.7)")
check(B.solve_asens({p: [0.7]}) == ODESLV.NORMAL and close(B.val_function_gradient(), [[-math.exp(-0.7)]], 1e-7),
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
C = ODESLV(G3)
C.add_domain(t, FFDom(0., 1., 4, FFDom.LGR, 3))
C.set_evolution_domain(t)
C.add_state(x, [t], ref=1.)
C.add_equation(OpP(x, t) + x, [t], [INT])
C.add_equation(x - 1., [t], [FFDom.LB], ii)
C.add_transition(x + 1., x, t, .5)
C.add_output(OpE(x, t, 1.))
C.options.RTOL = C.options.ATOL = 1e-10
check(C.setup(), "setup() with a transition  (%s)" % (C.extract_error() or "ok"))
check(C.solve([]) == ODESLV.NORMAL and close(C.val_function(), [math.exp(-1.) + math.exp(-0.5)], 1e-8),
      "x(1) == exp(-1) + exp(-0.5)")

# %% [markdown]
# ## D. A distributed control
#
# x' = u(t), x(0) = 0, u piecewise constant on 4 elements of [0, 1]:
# x(1) = sum of u_e / 4 and dx(1)/du_e = 1/4. The control values are given by a
# generator of the DOF index.

# %%
G4 = FFGraph()
t, x, u = G4.add_var("t"), G4.add_var("x(t)"), G4.add_var("u(t)")
D = ODESLV(G4)
D.add_domain(t, FFDom(0., 1., 4, FFDom.LGR, 3))
D.set_evolution_domain(t)
D.add_state(x, [t], ref=0.)
D.add_input(u, [t], type=FFDom.LG, n_node=1, ref=1.)
D.add_equation(OpP(x, t) - u, [t], [INT])
D.add_equation(x, [t], [FFDom.LB], ii)
D.add_output(OpE(x, t, 1.))
D.register_control(u)
D.options.RTOL = D.options.ATOL = D.options.RTOLS = D.options.ATOLS = 1e-10
check(D.setup(), "setup() with a distributed control  (%s)" % (D.extract_error() or "ok"))
check(D.n_control_dof() == 4, "4 control DOFs (one per element)")
gen = lambda dof: 1. + dof.element[t]                       # u_e = 1, 2, 3, 4
check(D.solve_fsens({u: gen}) == ODESLV.NORMAL and close(D.val_function(), [2.5], 1e-8),
      "x(1) == (1+2+3+4)/4 = 2.5  (generator of the DOF index)")
g = D.val_function_gradient(u)
print("  dx(1)/du =", g)
check(close(g, [[.25]] * 4, 1e-7), "dx(1)/du_e == 1/4 for every element (one row per DOF)")
P = D.set_parameter_values(u, [1., 2., 3., 4.], [0.] * D.np())
check(close(D.get_parameter_values(u, P), [1., 2., 3., 4.], 0.), "set/get_parameter_values round trip")

# %% [markdown]
# ## E. Options

# %%
o = ODESLV.Options()
print("  defaults: INTMETH", o.INTMETH, "| LINSOL", o.LINSOL, "| RTOL", o.RTOL, "| REDUCE.ORDER", o.REDUCE.ORDER)
check(isinstance(o, FFModel.Options), "ODESLV.Options extends FFModel.Options")
print("  doc of RTOL:", ODESLV.Options.RTOL.__doc__)

# %% [markdown]
# ## Summary

# %%
print("\ncheck_odeslv: %d passed, %d failed -- %s" % (npass, nfail, "FAILURES" if nfail else "ALL PASS"))
if __name__ == "__main__" and "ipykernel" not in sys.modules and nfail:
    sys.exit(1)
