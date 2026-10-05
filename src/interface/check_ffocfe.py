# %% [markdown]
# # CRONOS Python bindings -- `FFOCFESLV` and `FFOCFERES` check
#
# Exercises the two OCFESLV operations of the DAG: `FFOCFESLV` (a solve as a DAG
# function of its inputs, one- and two-map forms) against closed forms, its DAG
# identity, policies and options; and `FFOCFERES` (the collocation residuals and
# outputs as a DAG function of the collocation vectors) at the converged solution,
# with a finite-difference check of its Jacobian. Run as a script, it exits with
# status 1 on any failure.

# %%
import gc
import math
import os
import sys

for d in reversed(os.environ.get("CRONOS_PYPATH", "").split(":")):
    if d and d not in sys.path:
        sys.path.insert(0, d)

import numpy as np
from pymcpp import FFGraph, FFPartial, FFEval
from cronos import FFDom, FFModel, OCFESLV, FFOCFESLV, FFOCFERES

npass = nfail = 0


def check(cond, what):
    """Record and print one check."""
    global npass, nfail
    print("  %s  %s" % ("PASS" if cond else "FAIL", what))
    if cond:
        npass += 1
    else:
        nfail += 1


OpP, OpE = FFPartial(), FFEval()
Role = FFModel.EqnRole
P0, Q0 = 0.7, 1.3


def model(G, marching=False):
    """x' = -p q x, x(0) = 1 on t in [0,1]; output x(1) = exp(-p q)."""
    t, x, p, q = G.add_var("t"), G.add_var("x"), G.add_var("p"), G.add_var("q")
    S = OCFESLV(G)
    S.add_domain(t, FFDom(0., 1., 8, FFDom.LGR, 6))
    S.set_evolution_domain(t)
    S.add_state(x, [t], ref=1.)
    S.add_input(p, ref=P0)
    S.add_input(q, ref=Q0)
    S.add_equation(OpP(x, t) + p * q * x, [t], [FFDom.ALL - FFDom.LB], OCFESLV.EqnOptions(Role.INTERIOR))
    S.add_equation(x - 1., [t], [FFDom.LB], OCFESLV.EqnOptions(Role.INITIAL))
    S.add_output(OpE(x, t, 1.))
    S.options.INTERFACE.IMPOSITION = OCFESLV.Options.IC_STRONG
    S.options.SOLVE.MARCHING = marching
    S.options.SOLVE.RES_TOL = 1e-12
    S.options.DISPLAY_LEVEL = 0
    S.setup()
    return S, (t, x, p, q)


def jacobian(G, y, X, xv):
    rows, cols, dy = G.fdiff(y, X)
    vals = G.eval(dy, X, xv) if dy else []
    J = np.zeros((len(y), len(X)))
    for i, j, v in zip(rows, cols, vals):
        J[i, j] = v
    return J


E = math.exp(-P0 * Q0)

# %% [markdown]
# ## A. FFOCFESLV

# %%
for marching in (False, True):
    mode = "marching" if marching else "monolithic"
    G = FFGraph()
    S, (t, x, p, q) = model(G, marching)
    P, Q = G.add_var("P"), G.add_var("Q")
    Op = FFOCFESLV()
    # COPY: each embedding holds its own solver (and control registry), so both can coexist in one DAG
    y = Op({p: [P], q: [Q]}, S, FFOCFESLV.COPY)         # one-map: both inputs differentiated
    v = G.eval(y, [P, Q], [P0, Q0])
    check(len(y) == 1 and abs(v[0] - E) < 1e-9, "%s one-map: x(1) == exp(-pq)  (%.12f)" % (mode, v[0]))
    J = jacobian(G, y, [P, Q], [P0, Q0])
    check(np.allclose(J, [[-Q0 * E, -P0 * E]], atol=1e-8), "%s one-map: fdiff gradient == (-q, -p) exp(-pq)" % mode)
    y2 = Op({p: [P]}, {q: [Q]}, S, FFOCFESLV.COPY)      # two-map: q passes through
    check(str(y[0]) != str(y2[0]), "%s: one-map and two-map embeddings are separate operations (lt)" % mode)
    try:
        J2 = jacobian(G, y2, [P, Q], [P0, Q0])
        J1 = jacobian(G, y, [P, Q], [P0, Q0])
        check(np.allclose(J2, [[-Q0 * E, 0.]], atol=1e-8) and np.allclose(J1, [[-Q0 * E, -P0 * E]], atol=1e-8),
              "%s two-map: d/dQ == 0 (map 2), and the one-map gradient is unchanged" % mode)
    except Exception as e:
        check(False, "%s two-map derivatives (%s)" % (mode, str(e)[:70]))

# %% [markdown]
# ### Policies and options

# %%
G = FFGraph()
S, (t, x, p, q) = model(G)
P, Q = G.add_var("P"), G.add_var("Q")
y = FFOCFESLV()({p: [P], q: [Q]}, S, FFOCFESLV.SHALLOW)
del S
gc.collect()
check(abs(G.eval(y, [P, Q], [P0, Q0])[0] - E) < 1e-9, "SHALLOW: solver kept alive by the DAG after `del`")
# SHALLOW embeddings share the solver's control registry: a second embedding with other maps invalidates the first,
# which then refuses to evaluate (rather than return wrong derivatives) -- use COPY to embed one solver several times
G3 = FFGraph()
S3, (t3, x3, p3, q3) = model(G3)
P3, Q3 = G3.add_var("P"), G3.add_var("Q")
ya = FFOCFESLV()({p3: [P3], q3: [Q3]}, S3)
yb = FFOCFESLV()({p3: [P3]}, {q3: [Q3]}, S3)
try:
    jacobian(G3, ya, [P3, Q3], [P0, Q0])
    check(False, "SHALLOW: a second embedding with other maps invalidates the first (clear error)")
except RuntimeError as e:
    check("control registry changed" in str(e), "SHALLOW: a second embedding with other maps invalidates the first (clear error)")
S2, (t2, x2, p2, q2) = model(FFGraph())
try:
    FFOCFESLV()({p2: [P], q2: [Q]}, S2, FFOCFESLV.TRANSFER)
    check(False, "TRANSFER refused")
except Exception as e:
    check("TRANSFER" in str(e), "TRANSFER refused: %s" % str(e)[:60])
o = FFOCFESLV.options
o.NP2NF = 2.5
check(FFOCFESLV.options.NP2NF == 2.5 and FFOCFESLV.options.GRADIENT == FFOCFESLV.AUTO,
      "options are shared and writable (NP2NF, GRADIENT default AUTO)")
o.NP2NF = 3.

# %% [markdown]
# ### Repeated operands
#
# The same DAG variable fed to two degrees of freedom (an LGL input's shared element-end node, a heater
# parametrised in the DAG): the value, the derivative and the printed expression must match a distinct-operand
# embedding (MC++ FFOp::evaluate_external moved the repeated operand twice: an empty FFExpr, 2 Oct 2026).

# %%
G4 = FFGraph()
t4, x4, p4 = G4.add_var("t"), G4.add_var("x"), G4.add_var("p")
S4 = OCFESLV(G4)
S4.add_domain(t4, FFDom(0., 1., 2, FFDom.LGR, 4))
S4.set_evolution_domain(t4)
S4.add_state(x4, [t4], ref=1.)
S4.add_input(p4, [t4], ref=P0, type=FFDom.LGL, n_node=1)       # one DOF per time element
S4.add_equation(OpP(x4, t4) + p4 * x4, [t4], [FFDom.ALL - FFDom.LB], OCFESLV.EqnOptions(Role.INTERIOR))
S4.add_equation(x4 - 1., [t4], [FFDom.LB], OCFESLV.EqnOptions(Role.INITIAL))
S4.add_output(OpE(x4, t4, 1.))
S4.options.DISPLAY_LEVEL = 0
S4.setup()
A4, B4 = G4.add_var("A"), G4.add_var("B")
yr = FFOCFESLV()({p4: [A4, A4]}, S4, FFOCFESLV.COPY)               # the same variable twice
yd = FFOCFESLV()({p4: [A4, A4 + 0. * B4]}, S4, FFOCFESLV.COPY)     # distinct operands
vr, vd = G4.eval(yr, [A4], [P0])[0], G4.eval(yd, [A4, B4], [P0, 1.])[0]
rr, cr, dr = G4.fdiff(yr, [A4])
gr = G4.eval(dr, [A4], [P0])[0]
rd, cd, dd = G4.fdiff(yd, [A4])
gd = G4.eval(dd, [A4, B4], [P0, 1.])[0]
check(abs(vr - math.exp(-P0)) < 1e-4 and abs(vr - vd) < 1e-14 and abs(gr - gd) < 1e-12 and abs(gr + math.exp(-P0)) < 1e-4,
      "repeated operand: value and d/dA identical to distinct operands (%.2e, %.2e), ~ exp(-p) on this grid"
      % (abs(vr - vd), abs(gr - gd)))
check("( A, A )" in yr[0].str(), "repeated operand: the expression prints (%s)" % yr[0].str()[-14:])

# %% [markdown]
# ## B. FFOCFERES -- the collocation residuals and outputs as a DAG function

# %%
G = FFGraph()
S, (t, x, p, q) = model(G)
var, inp = S.init()
S.set_input_values(p, [P0], inp)
S.set_input_values(q, [Q0], inp)
S.solve(var, inp)
f_ref = S.val_functions()
X = [G.add_var("X%d" % i) for i in range(S.n_colloc_sta())]
U = [G.add_var("U%d" % i) for i in range(S.n_colloc_inp())]
R = FFOCFERES()(X, U, S)
ne, nf = S.n_colloc_eqn(), S.n_colloc_fct()
check(len(R) == ne + nf, "FFOCFERES: n_colloc_eqn() + n_colloc_fct() = %d outputs" % (ne + nf))
vals = np.array(G.eval(R, X + U, list(var) + list(inp)))
check(np.abs(vals[:ne]).max() < 1e-9, "residuals at the converged solution ~ 0 (max %.1e)" % np.abs(vals[:ne]).max())
check(abs(vals[ne] - f_ref[0]) < 1e-12, "outputs == val_functions() (%.12f)" % vals[ne])
# Jacobian column for one state coefficient, against finite differences
k, h = 3, 1e-7
rows, cols, dR = G.fdiff(R, X)
dv = G.eval(dR, X + U, list(var) + list(inp)) if dR else []
col = np.zeros(ne + nf)
for i, j, w in zip(rows, cols, dv):
    if j == k:
        col[i] = w
xh = list(var); xh[k] += h
fd = (np.array(G.eval(R, X + U, xh + list(inp))) - vals) / h
check(np.allclose(col, fd, atol=1e-5), "Jacobian column %d == finite difference (max diff %.1e)" % (k, np.abs(col - fd).max()))
try:
    FFOCFERES()(X[:-1], U, S)
    check(False, "a wrong-size state vector refused")
except Exception as e:
    check("expected" in str(e), "a wrong-size state vector refused")
Rs = FFOCFERES()(X, U, S, FFOCFERES.SHALLOW)
check(str(Rs[0]) != str(R[0]), "COPY and SHALLOW residual operations are separate (lt)")

# %% [markdown]
# ## Summary

# %%
print("\ncheck_ffocfe: %d passed, %d failed -- %s" % (npass, nfail, "FAILURES" if nfail else "ALL PASS"))
if __name__ == "__main__" and "ipykernel" not in sys.modules and nfail:
    sys.exit(1)
