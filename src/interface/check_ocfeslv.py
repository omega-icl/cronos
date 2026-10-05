# %% [markdown]
# # CRONOS Python bindings -- `OCFESLV` check
#
# Exercises `cronos.OCFESLV`: a PDE solved monolithic and marching against closed
# forms (the model of `MIXED_red`), the causality refusal, forward and adjoint
# sensitivities of an output with respect to a control, the input vectors
# (`set_input_values`, `sample_input`, `get_input_values`) and their guards, and the
# options. Run as a script, it exits with status 1 on any failure.

# %%
import math
import os
import sys

for d in reversed(os.environ.get("CRONOS_PYPATH", "").split(":")):
    if d and d not in sys.path:
        sys.path.insert(0, d)

import numpy as np
import pymcpp
from pymcpp import FFGraph, FFPartial, FFIntegral, FFEval
from cronos import FFDom, FFModel, OCFESLV

npass = nfail = 0


def check(cond, what):
    """Record and print one check."""
    global npass, nfail
    print("  %s  %s" % ("PASS" if cond else "FAIL", what))
    if cond:
        npass += 1
    else:
        nfail += 1


OpP, OpI, OpE = FFPartial(), FFIntegral(), FFEval()
Role = FFModel.EqnRole

# %% [markdown]
# ## A. Heat equation, monolithic and marching, against closed forms
#
# u_t = a u_zz on (0,1] x (0,1), u = 0 at z = 0 and z = 1, u(0,z) = sin(pi z):
# u = exp(-a pi^2 t) sin(pi z), a = 1/4. Outputs: u(1, 1/2) = exp(-a pi^2) and
# the double integral (2/pi)(1 - exp(-a pi^2)) / (a pi^2).

# %%
A = 0.25
E1 = math.exp(-A * math.pi**2)
II = 2. / math.pi * (1. - E1) / (A * math.pi**2)


def heat(marching):
    G = FFGraph()
    t, z, u = G.add_var("t"), G.add_var("z"), G.add_var("u(t,z)")
    S = OCFESLV(G)
    S.add_domain(t, FFDom([0., 0.4, 1.], FFDom.LGR, 8))
    S.add_domain(z, FFDom(0., 1., 4, FFDom.LGL, 7))
    S.set_evolution_domain(t)
    S.add_state(u, [t, z], ref=0.5)
    T_NO_LB, Z_INT = FFDom.ALL - FFDom.LB, FFDom.ALL - FFDom.LB - FFDom.UB
    S.add_equation(OpP(u, t) - A * OpP(u, {z: 2}), [t, z], [T_NO_LB, Z_INT], OCFESLV.EqnOptions(Role.INTERIOR))
    S.add_equation(u, [t, z], [T_NO_LB, FFDom.LB], OCFESLV.EqnOptions(Role.BOUNDARY))
    S.add_equation(u, [t, z], [T_NO_LB, FFDom.UB], OCFESLV.EqnOptions(Role.BOUNDARY))
    S.add_equation(u - pymcpp.sin(math.pi * z), [t, z], [FFDom.LB, FFDom.ALL], OCFESLV.EqnOptions(Role.INITIAL))
    S.add_output(OpE(u, {t: 1, z: 1}, {t: 1., z: .5}))
    S.add_output(OpI(u, {t: 1, z: 1}))
    S.options.SOLVE.MARCHING = marching
    S.options.SOLVE.RES_TOL = 1e-12
    S.options.DISPLAY_LEVEL = 0
    ok = S.setup()
    var, inp = S.init()
    rep = S.solve(var, inp)
    return S, ok, rep, var


for marching in (False, True):
    mode = "marching" if marching else "monolithic"
    S, ok, rep, var = heat(marching)
    f = S.val_functions()
    print("  %-10s %r  f = %s" % (mode, rep, f))
    check(ok and rep.converged and S.is_marching() == marching, "%s: setup, solve converged" % mode)
    check(abs(f[0] - E1) < 1e-6 and abs(f[1] - II) < 1e-6, "%s: u(1,1/2) and the double integral == closed forms (1e-6)" % mode)
check(isinstance(var, np.ndarray) and var.size == S.n_colloc_sta(), "init: var is a NumPy array of n_colloc_sta() entries")

# %% [markdown]
# ## B. The causality refusal
#
# The value of a STATE at a point of the evolution direction cannot feed an
# equation of the solve that produces it.

# %%
G = FFGraph()
t, x, s = G.add_var("t"), G.add_var("x(t)"), G.add_var("s")
B = OCFESLV(G)
B.add_domain(t, FFDom(0., 1., 2, FFDom.LGR, 4))
B.set_evolution_domain(t)
B.add_state(x, [t], ref=1.)
B.add_state(s, [], ref=.5)
B.add_equation(OpP(x, t) + .7 * x, [t], [FFDom.ALL - FFDom.LB], OCFESLV.EqnOptions(Role.INTERIOR))
B.add_equation(x - 1., [t], [FFDom.LB], OCFESLV.EqnOptions(Role.INITIAL))
B.add_equation(s - OpE(x, t, 1.), [], [], OCFESLV.EqnOptions(Role.INTERIOR))
B.add_output(s)
B.options.DISPLAY_LEVEL = 0
check(not B.setup(), "setup() refuses a state evaluation inside an equation  (%s)" % B.setup_status)

# %% [markdown]
# ## C. Sensitivities of an output with respect to a control
#
# x' = -p x, x(0) = 1, output x(1) = exp(-p); p a scalar control:
# dx(1)/dp = -exp(-p), by forward and adjoint sensitivity analysis.

# %%
G = FFGraph()
t, x, p = G.add_var("t"), G.add_var("x(t)"), G.add_var("p")
C = OCFESLV(G)
C.add_domain(t, FFDom(0., 1., 4, FFDom.LGR, 6))
C.set_evolution_domain(t)
C.add_state(x, [t], ref=1.)
C.add_input(p, ref=0.7)
C.add_equation(OpP(x, t) + p * x, [t], [FFDom.ALL - FFDom.LB], OCFESLV.EqnOptions(Role.INTERIOR))
C.add_equation(x - 1., [t], [FFDom.LB], OCFESLV.EqnOptions(Role.INITIAL))
C.add_output(OpE(x, t, 1.))
C.register_control(p)
C.options.SOLVE.RES_TOL = 1e-12
C.options.DISPLAY_LEVEL = 0
check(C.setup(), "setup() with a control")
var, inp = C.init()
C.set_input_values(p, [0.7], inp)
check(C.get_input_values(p, inp) == [0.7] and C.encode_controls(inp) == [0.7], "set / get_input_values, encode_controls")
for name, solve in (("forward", C.solve_fsens), ("adjoint", C.solve_asens)):
    v = var.copy()
    ok = solve(v, inp)
    J = C.sens_jacobian()
    print("  %s: sens_jacobian %s = %s" % (name, J.shape, J.ravel()))
    check(ok and J.shape == (1, 1) and abs(J[0, 0] + math.exp(-0.7)) < 1e-8, "%s: dx(1)/dp == -exp(-0.7) (1e-8)" % name)

# %% [markdown]
# ## D. Inputs: functions of the coordinates, and the array guards

# %%
G = FFGraph()
t, x, w = G.add_var("t"), G.add_var("x(t)"), G.add_var("w(t)")
D = OCFESLV(G)
D.add_domain(t, FFDom(0., 1., 4, FFDom.LGR, 3))
D.set_evolution_domain(t)
D.add_state(x, [t], ref=0.)
D.add_input(w, [t], ref=0.)
D.add_equation(OpP(x, t) - w, [t], [FFDom.ALL - FFDom.LB], OCFESLV.EqnOptions(Role.INTERIOR))
D.add_equation(x, [t], [FFDom.LB], OCFESLV.EqnOptions(Role.INITIAL))
D.add_output(OpE(x, t, 1.))
D.options.DISPLAY_LEVEL = 0
D.setup()
var, inp = D.init()
D.sample_input(w, lambda c: 2. * c[t], inp)          # w(t) = 2t  ->  x(1) = 1
rep = D.solve(var, inp)
check(rep.converged and abs(D.val_functions()[0] - 1.) < 1e-10, "sample_input w = 2t: x(1) == 1")
offsets = np.array(D.node_colloc(w)).ravel()          # node offsets within an element (the first one)
t_nodes = np.array([e * 0.25 + o for e in range(4) for o in offsets])   # element-major, node-inner
check(np.allclose(D.get_input_values(w, inp), 2. * t_nodes), "get_input_values == 2 t at the input's collocation nodes")
for what, call, err in [
        ("var of the wrong size", lambda: D.solve(np.zeros(D.n_colloc_sta() + 1), inp), ValueError),
        ("var of dtype int", lambda: D.solve(np.zeros(D.n_colloc_sta(), dtype=int), inp), TypeError),
        ("var as a list", lambda: D.solve([0.] * D.n_colloc_sta(), inp), TypeError),
        ("non-contiguous var", lambda: D.solve(np.zeros(2 * D.n_colloc_sta())[::2], inp), ValueError)]:
    try:
        call()
        check(False, "%s refused" % what)
    except err:
        check(True, "%s refused (%s)" % (what, err.__name__))

# %% [markdown]
# ## E. Options and equation options

# %%
o = OCFESLV.Options()
print("  SOLVE.MARCHING", o.SOLVE.MARCHING, "| SOLVE.RES_TOL", o.SOLVE.RES_TOL, "| REDUCE.ORDER", o.REDUCE.ORDER)
check(isinstance(o, FFModel.Options), "OCFESLV.Options extends FFModel.Options")
eo = OCFESLV.EqnOptions(Role.INTERFACE)
check(isinstance(eo, FFModel.EqnOptions) and not eo.participate_in_classification,
      "OCFESLV.EqnOptions extends FFModel.EqnOptions; an INTERFACE row never enters classification")

# %% [markdown]
# ## Summary

# %%
print("\ncheck_ocfeslv: %d passed, %d failed -- %s" % (npass, nfail, "FAILURES" if nfail else "ALL PASS"))
if __name__ == "__main__" and "ipykernel" not in sys.modules and nfail:
    sys.exit(1)
