# %% [markdown]
# # CRONOS Python bindings -- `FFODESLV` check
#
# Exercises `cronos.FFODESLV`: an `ODESLV` solve embedded in the DAG, its outputs
# evaluated (`FFGraph.eval`) and differentiated (`FFGraph.fdiff` / `bdiff`) with
# forward and adjoint sensitivities; the one-map and two-map forms; DOFs given as
# lists or by a generator; the copy policies; a model with a transition, numerically
# and symbolically (`SYMDIFF`). Checked against closed forms and the
# C++ driver `ODESLV_test1`. Run as a script, it exits with status 1 on any failure.

# %%
import gc
import math
import os
import sys

for d in reversed(os.environ.get("CRONOS_PYPATH", "").split(":")):
    if d and d not in sys.path:
        sys.path.insert(0, d)

from pymcpp import FFGraph, FFPartial, FFIntegral, FFEval
from cronos import FFDom, FFModel, ODESLV, FFODESLV

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
    if isinstance(a, (list, tuple)):
        return len(a) == len(b) and all(close(x, y, rtol) for x, y in zip(a, b))
    return abs(a - b) <= rtol * max(1., abs(b))


def jacobian(G, deps, indeps, vals, mode="fdiff"):
    """Dense Jacobian d deps / d indeps at vals, from FFGraph.fdiff / bdiff (sparse)."""
    rows, cols, dvars = getattr(G, mode)(deps, indeps)
    dv = G.eval(dvars, indeps, vals) if dvars else []
    J = [[0.] * len(indeps) for _ in deps]
    for r, c, v in zip(rows, cols, dv):
        J[r][c] = v
    return J


OpP, OpI, OpE = FFPartial(), FFIntegral(), FFEval()
INT = FFDom.ALL - FFDom.LB
ii = FFModel.EqnOptions(FFModel.EqnRole.INITIAL)

# %% [markdown]
# ## A. A closed form through the DAG
#
# x' = -p x, x(0) = 1: y = x(1) = exp(-P), dy/dP = -exp(-P), for both
# sensitivity modes and both DAG differentiation directions.

# %%
G = FFGraph()
t, x, p = G.add_var("t"), G.add_var("x(t)"), G.add_var("p")
B = ODESLV(G)
B.options.SOLVE.VERBOSE = 0
B.add_domain(t, FFDom(0., 1., 2, FFDom.LGR, 3))
B.set_evolution_domain(t)
B.add_state(x, [t], ref=1.)
B.add_input(p)
B.add_equation(OpP(x, t) + p * x, [t], [INT])
B.add_equation(x - 1., [t], [FFDom.LB], ii)
B.add_output(OpE(x, t, 1.))
for k in ("RTOL", "ATOL", "RTOLS", "ATOLS", "RTOLB", "ATOLB"):
    setattr(B.options.SOLVE, k, 1e-10)
check(B.setup(), "solver setup()")
P = G.add_var("P")
OpODE = FFODESLV()
y = OpODE({p: [P]}, B)                                # one-map form
check(len(y) == 1, "one output: x(1)  [%s]" % y[0])
check(close(G.eval(y, [P], [0.7]), [math.exp(-0.7)], 1e-8), "eval: x(1) == exp(-0.7)")
for mode in (FFODESLV.FORWARD, FFODESLV.ADJOINT):
    FFODESLV.options.GRADIENT = mode
    for d in ("fdiff", "bdiff"):
        J = jacobian(G, y, [P], [0.7], d)
        check(close(J, [[-math.exp(-0.7)]], 1e-7), "%s, %s: dx(1)/dP == -exp(-0.7)" % (mode.name, d))
FFODESLV.options.GRADIENT = FFODESLV.AUTO

# %% [markdown]
# ## B. The model of `ODESLV_test1`, against the C++ driver

# %%
G2 = FFGraph()
t = G2.add_var("t")
x0, x1 = G2.add_var("x0(t)"), G2.add_var("x1(t)")
p0, p1 = G2.add_var("p0"), G2.add_var("p1")
A = ODESLV(G2)
A.options.SOLVE.VERBOSE = 0
A.add_domain(t, FFDom(0., 10., 4, FFDom.LGR, 4))
A.add_state(x0, [t], ref=1.2)
A.add_state(x1, [t], ref=1.1)
A.add_input(p0)
A.add_input(p1)
A.set_evolution_domain(t)
A.add_equation(OpP(x0, t) - p0 * x0 * (1. - x1), [t], [INT])
A.add_equation(OpP(x1, t) - p0 * x1 * (x0 - 1.), [t], [INT])
A.add_equation(x0 - 1.2, [t], [FFDom.LB], ii)
A.add_equation(x1 - (1.1 + 0.01 * p1), [t], [FFDom.LB], ii)
A.add_output(OpI(x1, t))
A.add_output(OpE(x0 * x1, t, 10.))
o = A.options
o.SOLVE.INTMETH, o.SOLVE.NLINSOL, o.SOLVE.LINSOL = ODESLV.Options.MSBDF, ODESLV.Options.NEWTON, ODESLV.Options.SPARSE
o.SOLVE.NMAX = 2000
o.SOLVE.ATOL = o.SOLVE.ATOLB = o.SOLVE.ATOLS = o.SOLVE.RTOL = o.SOLVE.RTOLB = o.SOLVE.RTOLS = 1e-9
o.SOLVE.QERR = o.SOLVE.QERRS = True
o.SOLVE.ASACHKPT = 2000
check(A.setup(), "solver setup()")
P0, P1 = G2.add_var("P0"), G2.add_var("P1")
Y = OpODE({p0: [P0], p1: [P1]}, A)
F_REF = [10.05033, 0.8200303]
FP_REF = [[-0.7158895, 1.974535], [-0.001896783, 0.0007985682]]   # [direction][output]
check(close(G2.eval(Y, [P0, P1], [2.96, 3.]), F_REF, 2e-6), "eval == C++ outputs (7 digits)")
J = jacobian(G2, Y, [P0, P1], [2.96, 3.])
check(close(J, [list(c) for c in zip(*FP_REF)], 5e-6), "fdiff Jacobian == C++ gradient (transposed: [output][parameter])")

# %% [markdown]
# ## C. The two-map form
#
# Map 1 = {p0}: differentiated numerically; map 2 = {p1}: its value passes through,
# with a zero numerical derivative by contract. On a FRESH DAG: see the note below on
# a one-map and a two-map operation of the same solver over the same variables.

# %%
G2b = FFGraph()
t = G2b.add_var("t")
x0, x1 = G2b.add_var("x0(t)"), G2b.add_var("x1(t)")
p0, p1 = G2b.add_var("p0"), G2b.add_var("p1")
A = ODESLV(G2b)
A.options.SOLVE.VERBOSE = 0
A.add_domain(t, FFDom(0., 10., 4, FFDom.LGR, 4))
A.add_state(x0, [t], ref=1.2)
A.add_state(x1, [t], ref=1.1)
A.add_input(p0)
A.add_input(p1)
A.set_evolution_domain(t)
A.add_equation(OpP(x0, t) - p0 * x0 * (1. - x1), [t], [INT])
A.add_equation(OpP(x1, t) - p0 * x1 * (x0 - 1.), [t], [INT])
A.add_equation(x0 - 1.2, [t], [FFDom.LB], ii)
A.add_equation(x1 - (1.1 + 0.01 * p1), [t], [FFDom.LB], ii)
A.add_output(OpI(x1, t))
A.add_output(OpE(x0 * x1, t, 10.))
o = A.options
o.SOLVE.INTMETH, o.SOLVE.NLINSOL, o.SOLVE.LINSOL = ODESLV.Options.MSBDF, ODESLV.Options.NEWTON, ODESLV.Options.SPARSE
o.SOLVE.NMAX = 2000
o.SOLVE.ATOL = o.SOLVE.ATOLB = o.SOLVE.ATOLS = o.SOLVE.RTOL = o.SOLVE.RTOLB = o.SOLVE.RTOLS = 1e-9
o.SOLVE.QERR = o.SOLVE.QERRS = True
o.SOLVE.ASACHKPT = 2000
A.setup()
G2 = G2b
P0, P1 = G2.add_var("P0"), G2.add_var("P1")
Y2 = OpODE({p0: [P0]}, {p1: [P1]}, A)
check(close(G2.eval(Y2, [P0, P1], [2.96, 3.]), F_REF, 2e-6), "two-map: same outputs")
J2 = jacobian(G2, Y2, [P0, P1], [2.96, 3.])
check(close([r[0] for r in J2], [FP_REF[0][0], FP_REF[0][1]], 5e-6) and all(r[1] == 0. for r in J2),
      "two-map: d/dP0 == C++, d/dP1 == 0 (map 2)")
y1 = OpODE(1, {p0: [P0]}, {p1: [P1]}, A)
check(close(G2.eval([y1], [P0, P1], [2.96, 3.]), [F_REF[1]], 2e-6), "single-output form (idep = 1)")

# %% [markdown]
# A ONE-map operation of the same solver over the same variables is a DIFFERENT
# operation (FFBaseDE::lt compares the differentiated set), with its own
# derivatives; calls differing only by `idep` share one operation.

# %%
Y1 = OpODE({p0: [P0], p1: [P1]}, A)
check([str(v) for v in Y1] != [str(v) for v in Y2], "one-map after two-map: a separate operation")
check(close(jacobian(G2, Y1, [P0, P1], [2.96, 3.]), [list(c) for c in zip(*FP_REF)], 5e-6),
      "  ... with the full derivatives (d/dP1 nonzero)")
check(str(OpODE(0, {p0: [P0], p1: [P1]}, A)) == str(Y1[0]) and str(OpODE(1, {p0: [P0], p1: [P1]}, A)) == str(Y1[1]),
      "idep calls share the one-map operation")

# %% [markdown]
# ## D. A distributed control, its DOFs given by a generator
#
# x' = u(t), x(0) = 0, u piecewise constant on 4 elements: x(1) = sum(U_e) / 4.

# %%
G3 = FFGraph()
t, x, u = G3.add_var("t"), G3.add_var("x(t)"), G3.add_var("u(t)")
D = ODESLV(G3)
D.options.SOLVE.VERBOSE = 0
D.add_domain(t, FFDom(0., 1., 4, FFDom.LGR, 3))
D.set_evolution_domain(t)
D.add_state(x, [t], ref=0.)
D.add_input(u, [t], type=FFDom.LG, n_node=1, ref=1.)
D.add_equation(OpP(x, t) - u, [t], [INT])
D.add_equation(x, [t], [FFDom.LB], ii)
D.add_output(OpE(x, t, 1.))
for k in ("RTOL", "ATOL", "RTOLS", "ATOLS"):
    setattr(D.options.SOLVE, k, 1e-10)
check(D.setup(), "solver setup()")
U = [G3.add_var("U%d" % e) for e in range(4)]
yd = OpODE({u: [U[d.element[t]] for d in D.input_dofs(u)]}, D)   # by DOF index: element e -> U[e]
vals = [1., 2., 3., 4.]
check(close(G3.eval(yd, U, vals), [2.5], 1e-8), "x(1) == (1+2+3+4)/4  (variables by DOF index, a list over input_dofs)")
check(close(jacobian(G3, yd, U, vals), [[.25] * 4], 1e-7), "dx(1)/dU_e == 1/4")
amp = G3.add_var("amp")                                        # (not A: that name is the solver of sections B, E)
ya = OpODE({u: lambda coord: amp * coord[t]}, D)               # a profile parameterised by amp: u = amp t at the midpoints
check(close(G3.eval(ya, [amp], [2.]), [1.], 1e-8), "a callable is a function of the coordinates: u = amp t, x(1) == amp/2")
check(close(jacobian(G3, ya, [amp], [2.]), [[.5]], 1e-7), "dx(1)/damp == 1/2")

# %% [markdown]
# ## E. Copy policies
#
# SHALLOW refers to the solver, which the binding keeps alive as long as the DAG;
# TRANSFER would hand Python's object to the DAG and is refused.

# %%
G4 = FFGraph()
t, x, p = G4.add_var("t"), G4.add_var("x(t)"), G4.add_var("p")
S4 = ODESLV(G4)
S4.options.SOLVE.VERBOSE = 0
S4.add_domain(t, FFDom(0., 1., 2, FFDom.LGR, 3))
S4.set_evolution_domain(t)
S4.add_state(x, [t], ref=1.)
S4.add_input(p)
S4.add_equation(OpP(x, t) + p * x, [t], [INT])
S4.add_equation(x - 1., [t], [FFDom.LB], ii)
S4.add_output(OpE(x, t, 1.))
S4.options.SOLVE.RTOL = S4.options.SOLVE.ATOL = 1e-10
S4.setup()
P4 = G4.add_var("P")
ys = OpODE({p: [P4]}, S4, FFODESLV.SHALLOW)
del S4
gc.collect()
check(close(G4.eval(ys, [P4], [0.3]), [math.exp(-0.3)], 1e-8), "SHALLOW: solver kept alive by the DAG after `del`")
try:
    OpODE({p0: [P0], p1: [P1]}, A, FFODESLV.TRANSFER)
    check(False, "TRANSFER refused")
except ValueError as e:
    check(True, "TRANSFER refused: %s" % str(e)[:50])

# %% [markdown]
# ## F. A transition, differentiated numerically and symbolically
#
# An impulsive dose: x' = -k x, x(0) = x0, x(0.4+) = x(0.4-) + d, so that
# x(1) = (x0 exp(-0.4 k) + d) exp(-0.6 k). The numeric route takes the
# sensitivities across the jump; the symbolic route (`SYMDIFF`) differentiates
# the MODEL -- the product carries the transition and, per direction, its
# derivative s(0.4+) = s(0.4-) + d_d -- and embeds the product in the DAG.

# %%
G5 = FFGraph()
t, x, k, d, x0 = (G5.add_var(n) for n in ("t", "x(t)", "k", "d", "x0"))
T5 = ODESLV(G5)
T5.options.SOLVE.VERBOSE = 0
T5.add_domain(t, FFDom(0., 1., 2, FFDom.LGR, 4))
T5.set_evolution_domain(t)
T5.add_state(x, [t], ref=1.)
for v in (k, d, x0):
    T5.add_input(v)
T5.add_equation(OpP(x, t) + k * x, [t], [INT])
T5.add_equation(x - x0, [t], [FFDom.LB], ii)
T5.add_transition(x + d, x, t, .4)
T5.add_output(OpE(x, t, 1.))
for o in ("RTOL", "ATOL", "RTOLS", "ATOLS", "RTOLB", "ATOLB"):
    setattr(T5.options.SOLVE, o, 1e-11)
check(T5.setup(), "setup() with a transition  (%s)" % (T5.extract_error() or "ok"))
V5 = [G5.add_var(n) for n in ("K", "D", "X0")]
KV, DV, XV = 1.3, 0.7, 1.1
y5 = OpODE({k: [V5[0]], d: [V5[1]], x0: [V5[2]]}, T5)
Y_REF = (XV * math.exp(-.4 * KV) + DV) * math.exp(-.6 * KV)
J_REF = [[-XV * math.exp(-KV) - .6 * DV * math.exp(-.6 * KV), math.exp(-.6 * KV), math.exp(-KV)]]
check(close(G5.eval(y5, V5, [KV, DV, XV]), [Y_REF], 1e-8), "eval: x(1) == closed form, across the jump")
check(close(jacobian(G5, y5, V5, [KV, DV, XV]), J_REF, 1e-7), "numeric (SYMDIFF empty): dx(1)/d(k, d, x0) == closed form")
FFODESLV.options.SYMDIFF = V5
try:
    for dd in ("fdiff", "bdiff"):
        check(close(jacobian(G5, y5, V5, [KV, DV, XV], dd), J_REF, 1e-7),
              "symbolic (SYMDIFF = all), %s: dx(1)/d(k, d, x0) == closed form" % dd)
except RuntimeError as e:
    check(False, "symbolic (SYMDIFF = all): %s" % str(e)[:80])
FFODESLV.options.SYMDIFF = []

# %% [markdown]
# ## G. Options (class-level)

# %%
print("  FFODESLV.options: GRADIENT", FFODESLV.options.GRADIENT, "| NP2NF", FFODESLV.options.NP2NF)
FFODESLV.options.NP2NF = 5.
check(FFODESLV.options.NP2NF == 5., "options are shared and writable")
FFODESLV.options.NP2NF = 3.

# %% [markdown]
# ## Summary

# %%
print("\ncheck_ffodeslv: %d passed, %d failed -- %s" % (npass, nfail, "FAILURES" if nfail else "ALL PASS"))
if __name__ == "__main__" and "ipykernel" not in sys.modules and nfail:
    sys.exit(1)
