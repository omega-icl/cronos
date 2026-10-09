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
# ## add_output index and blk_fct

# %%
_FA = FFGraph(); _t7, _z7, _u7 = _FA.add_var("t"), _FA.add_var("z"), _FA.add_var("u")
_S7 = OCFESLV(_FA)
_S7.add_domain(_t7, FFDom(0., 1., 2, FFDom.LGR, 3)); _S7.add_domain(_z7, FFDom(0., 1., 2, FFDom.LGL, 4)); _S7.set_evolution_domain(_t7)
_S7.add_state(_u7, [_t7, _z7], ref=0.)
_NLB = FFDom.ALL - FFDom.LB; _INN = FFDom.ALL - FFDom.LB - FFDom.UB; _P = FFPartial(); _R = OCFESLV.EqnRole
_S7.add_equation(_P(_u7, _t7) - 0.1 * _P(_u7, {_z7: 2}), [_t7, _z7], [_NLB, _INN], OCFESLV.EqnOptions(_R.INTERIOR))
_S7.add_equation(_u7, [_t7, _z7], [_NLB, FFDom.LB], OCFESLV.EqnOptions(_R.BOUNDARY))
_S7.add_equation(_u7 - 1., [_t7, _z7], [_NLB, FFDom.UB], OCFESLV.EqnOptions(_R.BOUNDARY))
_S7.add_equation(_u7, [_t7, _z7], [FFDom.LB, FFDom.ALL], OCFESLV.EqnOptions(_R.INITIAL))
_j0 = _S7.add_output(_u7, [_t7, _z7], point=[1., .5])
_j1 = _S7.add_output(_u7, [_z7], masks=[FFDom.ALL], at={_t7: 1.})
_j2 = _S7.add_output(_u7, [_t7, _z7], point=[.5, .5])
_S7.options.DISPLAY_LEVEL = 0
check((_j0, _j1, _j2) == (0, 1, 2), "add_output returns the index: 0, 1, 2 (%s)" % str((_j0, _j1, _j2)))
check(_S7.setup(), "the model with two point outputs and a profile sets up")
_bb = [_S7.blk_fct(k) for k in (_j0, _j1, _j2)]
check(_bb[0] == (0, 1) and _bb[1][0] == 1 and _bb[1][1] > 1 and _bb[2] == (1 + _bb[1][1], 1),
      "blk_fct(index): the point outputs take one value, the profile its nodes, consecutively (%s)" % str(_bb))

# %% [markdown]
# ## AUTO.HYP_CLOSURE off: an open outflow face is refused (2026-10-07)
#
# Scalar advection u_t + a u_z = 0 (a = 1), exact u = (z - a t)^2: inflow condition at z = 0, outflow at z = 1.  With the
# automatic closure on, the outflow face is closed and the solution exact; with it off and the face left open, the
# system would be under-determined -- setup() refuses (status HYP_CLOSURE_MISSING) instead of solving it silently wrong.

# %%
def advection(hyp_closure):
    G = FFGraph(); t, z, u, a = G.add_var("t"), G.add_var("z"), G.add_var("u"), G.add_var("a")
    P, E, R = FFPartial(), FFEval(), OCFESLV.EqnRole
    NLB, INN = FFDom.ALL - FFDom.LB, FFDom.ALL - FFDom.LB - FFDom.UB
    S = OCFESLV(G)
    S.add_domain(t, FFDom(0., .25, 2, FFDom.LGR, 3)); S.add_domain(z, FFDom(0., 1., 4, FFDom.LGL, 4))
    S.set_evolution_domain(t)
    S.add_state(u, [t, z], ref=0.); S.add_input(a, ref=1.)
    S.add_equation(P(u, t) + a * P(u, z), [t, z], [NLB, INN], OCFESLV.EqnOptions(R.INTERIOR))
    S.add_equation(u - z * z, [t, z], [FFDom.LB, FFDom.ALL], OCFESLV.EqnOptions(R.INITIAL))
    S.add_equation(u - a * a * t * t, [t, z], [NLB, FFDom.LB], OCFESLV.EqnOptions(R.BOUNDARY))   # inflow at z = 0
    S.add_output(E(u, {t: 1, z: 1}, {t: .25, z: .5}))
    S.options.AUTO.HYP_CLOSURE = hyp_closure
    S.options.SOLVE.MARCHING = False; S.options.DISPLAY_LEVEL = 0
    return S

_H1 = advection(True)
_ok1 = _H1.setup()
_v1 = float("nan")
if _ok1:
    _var, _inp = _H1.init(); _rep = _H1.solve(_var, _inp); _v1 = _H1.val_functions()[0]
check(_ok1 and abs(_v1 - .25 ** 2) < 1e-8,
      "AUTO.HYP_CLOSURE on: the outflow face is closed, u(.25,.5) exact (|err| = %.1e)" % abs(_v1 - .25 ** 2))
_H0 = advection(False)
_ok0 = _H0.setup()
check(not _ok0 and _H0.setup_status.name == "HYP_CLOSURE_MISSING",
      "AUTO.HYP_CLOSURE off, outflow face open: setup() refuses (%s)" % _H0.setup_status)

# %% [markdown]
# ### Backends requested explicitly
# `SOLVE_SPQR`, `DET_SPQR` and `DET_EIGEN` exist in every build (2026-10-09).  Requested in a build without that
# backend, setup() refuses with `BACKEND_UNAVAILABLE` -- it used to be a compile error (SOLVE_SPQR) or a silently
# skipped check (DET_*).  In a build with it, the same model sets up as usual.

# %%
import cronos as _cr
_be = _cr.build_info()["backends"]
check(hasattr(OCFESLV.Options, "SOLVE_SPQR") and hasattr(OCFESLV.Options, "DET_SPQR")
      and hasattr(OCFESLV.Options, "DET_EIGEN"), "SOLVE_SPQR, DET_SPQR, DET_EIGEN exist in every build")
for _what, _field, _value, _have in (
        ("SOLVE.FACTORIZATION = SOLVE_SPQR", "SOLVE", OCFESLV.Options.SOLVE_SPQR, _be["SPQR"]),
        ("DETERMINACY.BACKEND = DET_SPQR", "DETERMINACY", OCFESLV.Options.DET_SPQR, _be["SPQR"]),
        ("DETERMINACY.BACKEND = DET_EIGEN", "DETERMINACY", OCFESLV.Options.DET_EIGEN, _be["Eigen"])):
    _B = advection(True)
    if _field == "SOLVE": _B.options.SOLVE.FACTORIZATION = _value
    else:                 _B.options.DETERMINACY.BACKEND = _value
    _okb = _B.setup()
    if _have:
        check(_okb and _B.setup_status.name == "OK", "%s, backend built in: setup() succeeds (%s)" % (_what, _B.setup_status))
    else:
        check(not _okb and _B.setup_status.name == "BACKEND_UNAVAILABLE",
              "%s, backend NOT built in: setup() refuses (%s)" % (_what, _B.setup_status))

# %% [markdown]
# ## Summary

# %%
print("\ncheck_ocfeslv: %d passed, %d failed -- %s" % (npass, nfail, "FAILURES" if nfail else "ALL PASS"))
if __name__ == "__main__" and "ipykernel" not in sys.modules and nfail:
    sys.exit(1)
