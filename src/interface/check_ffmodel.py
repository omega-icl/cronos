# %% [markdown]
# # CRONOS Python bindings -- `FFModel` check
#
# Exercises the `cronos.FFModel` binder: model declaration, `setup()`, the setup
# report and the model report, the declared records (`t_Eqn`, `t_Fct`), options,
# inputs and controls, transitions, and the refusals. Every check prints `PASS` or
# `FAIL`; run as a script (`python check_ffmodel.py`) it exits with status 1 on any
# failure, so it can also serve as a gate.
#
# Requires the `pymcpp` module (with `mc_ocbase`) and the `cronos` module on the
# Python path -- set `CRONOS_PYPATH` (colon-separated directories) or edit the cell
# below.

# %%
import os
import sys

for d in reversed(os.environ.get("CRONOS_PYPATH", "").split(":")):
    if d and d not in sys.path:
        sys.path.insert(0, d)

import pymcpp
import cronos
from pymcpp import FFGraph, FFMon, FFPartial, FFIntegral, FFEval
from cronos import FFDom, FFModel

print("pymcpp:", pymcpp.__file__)
print("cronos:", cronos.__file__)

npass = nfail = 0


def check(cond, what):
    """Record and print one check."""
    global npass, nfail
    print("  %s  %s" % ("PASS" if cond else "FAIL", what))
    if cond:
        npass += 1
    else:
        nfail += 1


def names(d):
    """A dict keyed by FFVar, re-keyed by variable name (FFVar hashes by identity)."""
    return {str(k): v for k, v in d.items()}


OpP, OpI, OpE = FFPartial(), FFIntegral(), FFEval()
INT = FFDom.ALL - FFDom.LB  # (a,b]: the evolution interior with the end point

# %% [markdown]
# ## A. An ODE model, its setup report and its model report
#
# The model of `ODESLV_test1`: two states, two parameters, an integral output and a
# terminal point output, on 4 Radau elements.

# %%
G = FFGraph()
t = G.add_var("t")
x0, x1 = G.add_var("x0(t)"), G.add_var("x1(t)")
p0, p1 = G.add_var("p0"), G.add_var("p1")

A = FFModel(G)
A.add_domain(t, FFDom(0., 10., 4, FFDom.LGR, 4))
A.set_evolution_domain(t)
A.add_state(x0, [t], ref=1.2)
A.add_state(x1, [t], ref=1.1)
A.add_input(p0)
A.add_input(p1)
io = FFModel.EqnOptions(FFModel.EqnRole.INTERIOR)
ii = FFModel.EqnOptions(FFModel.EqnRole.INITIAL)
A.add_equation(OpP(x0, t) - p0 * x0 * (1. - x1), [t], [INT], io)
A.add_equation(OpP(x1, t) - p0 * x1 * (x0 - 1.), [t], [INT], io)
A.add_equation(x0 - 1.2, [t], [FFDom.LB], ii)
A.add_equation(x1 - (1.1 + 0.01 * p1), [t], [FFDom.LB], ii)
A.add_output(OpI(x1, t))
A.add_output(OpE(x0 * x1, t, 10.))

print("declared equations (t_Eqn):")
for e in A.var_equation:
    print("   %-9s %s" % (e.opt.role.name, e))
print("declared outputs (t_Fct):")
for f in A.var_output:
    print("   %-9s %s" % (f.kind.name, f))

check([e.opt.role for e in A.var_equation] == [FFModel.EqnRole.INTERIOR] * 2 + [FFModel.EqnRole.INITIAL] * 2,
      "equation roles as declared")
check([names(e.dom) for e in A.var_equation] == [{"t": INT}] * 2 + [{"t": FFDom.LB}] * 2,
      "equation masks as declared ((a,b] and a)")
check(all(f.kind == FFModel.FctKind.POINT for f in A.var_output), "two scalar outputs")
# Before setup() the records hold the user's own variables: dicts returned by C++ can be indexed with them
# (FFVar compares and hashes by DAG identity). After setup() they hold the WORKING model's copies (another DAG).
k = next(iter(A.var_equation[0].dom))
check(k == t and hash(k) == hash(t) and A.var_equation[0].dom[t] == INT and A.var_equation[2].dom[t] == FFDom.LB,
      "returned dicts are indexed with the user's own FFVar (FFVar __eq__ / __hash__)")

# %% [markdown]
# Options: every field is bound, with the C++ documentation as its docstring.
# `DISPLAY_LEVEL = 1` makes `setup()` print its report; the binder redirects C++'s
# `std::cout` / `std::cerr` to Python, so the report appears here.

# %%
o = A.options
print("DISPLAY_LEVEL", o.DISPLAY_LEVEL, "| REDUCE.ORDER", o.REDUCE.ORDER,
      "| CLASSIFY.MODE", o.CLASSIFY.MODE, "| TTOL", o.TTOL)
print("doc of REDUCE.ORDER:", type(o.REDUCE).ORDER.__doc__)
A.options.DISPLAY_LEVEL = 1
ok = A.setup()
check(ok and A.setup_status == FFModel.SetupStatus.OK, "setup() succeeds  (%s)" % A.setup_status)
check(A.is_setup(), "is_setup() after setup")

# %%
A.report()
rep = str(A)
check("DIFFERENTIAL_ORDINARY" in rep and "AS AN ODE/DAE" not in rep, "report: classified DIFFERENTIAL_ORDINARY (no ODE/DAE preview section)")
check("(balanced)" in rep, "report: degrees of freedom balanced")
check("DEFERRED VALUES" not in rep and "holds a deferred value" not in rep and "OUTPUTS (2)" in rep,
      "report: the two outputs as declared, without the deferred-value plumbing")
# every declared state is listed -- also one a deferred value is taken of (2026-10-05: a point output of a state at
# the end of the evolution direction hid that state from the report)
_st = rep[rep.index("STATES"):rep.index("\n\n", rep.index("STATES"))]
check(all(("\n  %s " % str(v)) in _st for v in A.var_state),
      "report: every state of the model is listed (%d)" % len(A.var_state))
check(len(A.var_state) == 2, "two states")

# %% [markdown]
# ## B. A PDE with every output form
#
# The heat equation u_t = a u_zz on (0,1] x (0,1), u = 0 at z = 0 and z = 1,
# u(0,z) = 4 z (1-z) -- after `MIXED_red`, with a polynomial initial condition (see the
# note on `pymcpp.sin` in the README) -- with a point output, a one-sided
# point output, a distributed output, a distributed output at a fixed coordinate,
# and a double integral.

# %%
import math

G2 = FFGraph()
t, z, u = G2.add_var("t"), G2.add_var("z"), G2.add_var("u(t,z)")
B = FFModel(G2)
B.add_domain(t, FFDom([0., 0.4, 1.], FFDom.LGR, 8))
B.add_domain(z, FFDom(0., 1., 4, FFDom.LGL, 7))
B.set_evolution_domain(t)
B.add_state(u, [t, z], ref=lambda c: 4. * c[z] * (1. - c[z]))  # a reference FUNCTION of the coordinates
Z_INT = FFDom.ALL - FFDom.LB - FFDom.UB
B.add_equation(OpP(u, t) - 0.25 * OpP(u, {z: 2}), [t, z], [INT, Z_INT],
               FFModel.EqnOptions(FFModel.EqnRole.INTERIOR))
B.add_equation(u, [t, z], [INT, FFDom.LB], FFModel.EqnOptions(FFModel.EqnRole.BOUNDARY))
B.add_equation(u, [t, z], [INT, FFDom.UB], FFModel.EqnOptions(FFModel.EqnRole.BOUNDARY))
B.add_equation(u - 4. * z * (1. - z), [t, z], [FFDom.LB, FFDom.ALL],
               FFModel.EqnOptions(FFModel.EqnRole.INITIAL))
B.add_output(OpE(u, {t: 1, z: 1}, {t: 1., z: .5}))                # point value (operator form)
B.add_output(u, [t, z], point=[.4, .5], side=[FFDom.PLUS, FFDom.MINUS])  # one-sided point output
B.add_output(u, [t, z], masks=[FFDom.ALL, FFDom.ALL])            # distributed over both domains
B.add_output(u, [z], masks=[FFDom.ALL], at={t: 1.})              # distributed in z at t = 1
B.add_output(OpI(u, {t: 1, z: 1}))                                # double integral

kinds = [f.kind for f in B.var_output]
for f in B.var_output:
    print("   %-11s %s" % (f.kind.name, f))
check(kinds == [FFModel.FctKind.POINT] * 2 + [FFModel.FctKind.DISTRIBUTED] * 2 + [FFModel.FctKind.POINT],
      "output kinds: 2 point, 2 distributed, 1 integral (scalar)")
f1 = B.var_output[1]
check(names(f1.point) == {"t": .4, "z": .5} and names(f1.side) == {"t": 1}, "one-sided point output: t=0.4^+, z=0.5")
check(names(B.var_output[3].point) == {"t": 1.} and names(B.var_output[3].grid) == {"z": FFDom.ALL},
      "distributed output at t = 1: grid z, point t")
B.options.DISPLAY_LEVEL = 1
check(B.setup(), "setup() succeeds  (%s)" % B.setup_status)
B.report()

# %% [markdown]
# ## C. Inputs and controls
#
# A distributed input p(t) with a reference function, registered as a control: its
# degrees of freedom are the collocation values, one per node and element, located
# by `control_dofs` as `DofIndex` records. A second input w(t) is held at known
# values with `fix_input`.

# %%
G3 = FFGraph()
t, x, p, w = G3.add_var("t"), G3.add_var("x(t)"), G3.add_var("p(t)"), G3.add_var("w(t)")
C = FFModel(G3)
C.add_domain(t, FFDom(0., 1., 3, FFDom.LGR, 2))
C.set_evolution_domain(t)
C.add_state(x, [t], ref=1.)
C.add_input(p, [t], ref=lambda c: 1. + c[t])          # reference function of the coordinate
C.add_input(w, [t], continuity=[FFModel.CONTINUOUS_C0], ref=.5)
C.add_equation(OpP(x, t) + x - p - w, [t], [INT])
C.add_equation(x - 1., [t], [FFDom.LB])
C.add_output(OpE(x, t, 1.))
C.register_control(p)
check(C.setup(), "setup() succeeds  (%s)" % C.setup_status)
n = C.n_control_dof()
check(n == 3 * 2, "control DOFs of p(t): 3 elements x 2 nodes = %d" % n)
dofs = C.control_dofs(p)
# dicts returned by C++ are keyed by NEW FFVar objects; FFVar compares and hashes by DAG identity, so
# they can be indexed with the user's own variables
print("DOF layout of p(t) (element, node):", [(d.element[t], d.node[t]) for d in dofs])
check([(d.element[t], d.node[t]) for d in dofs] == [(e, k) for e in range(3) for k in range(2)],
      "DOFs ordered element-major, node-minor")
roles = {str(e.var): e.opt.role for e in C.var_equation}
check(list(roles.values()) == [FFModel.EqnRole.INTERIOR, FFModel.EqnRole.INITIAL]
      and all(e.opt.role_auto for e in C.var_equation),
      "AUTO roles: interior (a,b] -> INTERIOR; LB alone of the EVOLUTION domain -> INITIAL (role_auto recorded)")
check(B.var_equation[1].opt.role == FFModel.EqnRole.BOUNDARY,
      "an explicit BOUNDARY row on a SPATIAL bound keeps its role")

# %% [markdown]
# ## D. A transition
#
# x' = -x on (0,1], x(0) = 1, with a jump x(0.5^+) = x(0.5^-) + 1, declared as
# `left(tau^-) = right(tau^+)`. A bare `FFModel` only analyses the model: it validates
# the transition and lists it in the report (`ODESLV` and `OCFESLV` integrate it -- see
# their checks).

# %%
import contextlib
import io

G4 = FFGraph()
t, x = G4.add_var("t"), G4.add_var("x(t)")
D = FFModel(G4)
D.add_domain(t, FFDom(0., 1., 4, FFDom.LGR, 3))
D.set_evolution_domain(t)
D.add_state(x, [t], ref=1.)
D.add_equation(OpP(x, t) + x, [t], [INT])
D.add_equation(x - 1., [t], [FFDom.LB])
D.add_transition(x + 1., x, t, .5)      # x(0.5^-) + 1 = x(0.5^+)
D.add_output(OpE(x, t, 1.))
D.options.DISPLAY_LEVEL = 1
buf = io.StringIO()
with contextlib.redirect_stderr(buf), contextlib.redirect_stdout(buf):
    ok = D.setup()
msg = buf.getvalue()
check(ok and D.setup_status == FFModel.SetupStatus.OK, "setup() with a transition succeeds  (%s)" % D.setup_status)
rep = str(D)
print("\n".join(l for l in rep.split("\n")[rep.split("\n").index(next(l for l in rep.split("\n") if l.startswith("TRANSITIONS"))):][:2]))
check("TRANSITIONS (1)" in rep and "x(t) + 1 = x(t)" in rep and "validated at setup" in rep,
      "report lists the validated transition")

# %% [markdown]
# ## E. Refusals and errors
#
# * an evolution-direction evaluation of a STATE inside an equation is refused at
#   setup (causality) -- the refusal message comes through the redirected std::cerr;
# * Python-side keyword misuse raises `ValueError`;
# * C++ exceptions surface as `RuntimeError` (for `FFModel`) and as the exception
#   class of `FFDom`.

# %%
import contextlib
import io

G5 = FFGraph()
t, x, s = G5.add_var("t"), G5.add_var("x(t)"), G5.add_var("s")
E = FFModel(G5)
E.add_domain(t, FFDom(0., 1., 2, FFDom.LGR, 3))
E.set_evolution_domain(t)
E.add_state(x, [t], ref=1.)
E.add_state(s, [], ref=.5)
E.add_equation(OpP(x, t) + .7 * x, [t], [INT])
E.add_equation(x - 1., [t], [FFDom.LB])
E.add_equation(s - OpE(x, t, 1.))         # refused: a state's value at t = 1 inside an equation
buf = io.StringIO()
with contextlib.redirect_stderr(buf), contextlib.redirect_stdout(buf):
    ok = E.setup()
msg = buf.getvalue()
print(msg.strip()[:400])
check(not ok and E.setup_status != FFModel.SetupStatus.OK, "setup() refuses the model  (%s)" % E.setup_status)
check("EQUATION REFUSED" in msg and "STATE" in msg, "refusal message captured from C++ (redirected std::cerr)")

for what, call in [
        ("add_output: point and masks together", lambda: E.add_output(x, [t], point=[.5], masks=[0])),
        ("add_output: side without a point output", lambda: E.add_output(x, [t], masks=[0], side=[1])),
        ("add_input: n_node without type", lambda: E.add_input(G5.add_var("q"), [t], n_node=3)),
        ("add_input: continuity with a callable ref", lambda: E.add_input(G5.add_var("r"), [t], continuity=[FFModel.CONTINUOUS_C0], ref=lambda c: 1.))]:
    try:
        call()
        check(False, "ValueError for %s" % what)
    except ValueError:
        check(True, "ValueError for %s" % what)

try:
    E.add_output(x, [t], point=[.5, .2])   # 2 coordinates for 1 domain: C++ Exceptions::INDEX
    check(False, "C++ INDEX error raised")
except RuntimeError as e:
    check(True, "C++ INDEX error -> RuntimeError: %s" % str(e)[:60])

try:
    FFDom(1., 0., 2)                       # upper bound below lower bound
    check(False, "FFDom refuses inverted bounds")
except Exception as e:
    check(True, "FFDom refuses inverted bounds -> %s" % type(e).__name__)

# %% [markdown]
# ## add_output returns the output's index

# %%
from pymcpp import FFGraph as _G, FFEval as _E
_G5 = _G(); _t5, _z5, _u5 = _G5.add_var("t"), _G5.add_var("z"), _G5.add_var("u")
_M5 = FFModel(_G5)
_M5.add_domain(_t5, FFDom(0., 1., 2)); _M5.add_domain(_z5, FFDom(0., 1., 2, FFDom.LGL, 4))
_M5.add_state(_u5, [_t5, _z5])
_i0 = _M5.add_output(_E()(_u5, {_t5: 1, _z5: 1}, {_t5: 1., _z5: .5}))      # a scalar output
_i1 = _M5.add_output(_u5, [_t5, _z5], masks=[FFDom.ALL, FFDom.UB])         # a distributed output
_i2 = _M5.add_output(_u5, [_t5, _z5], point=[1., .5])                      # a point output
check((_i0, _i1, _i2) == (0, 1, 2), "add_output returns the output's index: declaration order from 0 (%s)" % str((_i0, _i1, _i2)))
check(len(_M5.var_output) == 3 and all(isinstance(i, int) for i in (_i0, _i1, _i2)), "the indices are ints, one per output")

# %% [markdown]
# ## AUTO.HYP_CLOSURE: an option again (2026-10-07)

# %%
_M9 = FFModel(FFGraph())
check(_M9.options.AUTO.HYP_CLOSURE is True, "AUTO.HYP_CLOSURE exists and defaults on")
check(hasattr(FFModel.SetupStatus, "HYP_CLOSURE_MISSING"),
      "setup status HYP_CLOSURE_MISSING exists (off with an outflow face left open is refused)")

# %% [markdown]
# ## Summary

# %%
print("\ncheck_ffmodel: %d passed, %d failed -- %s" % (npass, nfail, "FAILURES" if nfail else "ALL PASS"))
if __name__ == "__main__" and "ipykernel" not in sys.modules and nfail:
    sys.exit(1)
