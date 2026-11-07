#!/usr/bin/env python3
"""check_daeslv_mol.py -- a model distributed in space on DAESLV, from Python (MoL phases 7'b-6c, 7).

The declared model is the PDE; DAESLV integrates its semi-discretisation by the method of lines. The checks
mirror the C++ gates (DAESLV_molpde, DAESLV_moltime) with closed forms where one exists:
  (A) a rod with a Dirichlet inlet p: DAESLV's output = OCFESLV's, forward = adjoint = central difference in p;
  (B) a TIME-VARYING inlet g(t) collocated in t: one gradient row per declared DOF, = central differences;
  (C) a TRANSITION u(tau+) = u(tau-)/2 on a zero-flux rod with a source: J and dJ/dp = closed forms, exact on
      any mesh (a uniform field is in the operator's kernel);
  (D) the semi-discrete system's unknowns and the fields' nodes, from semidiscrete_states / semidiscrete_nodes.
"""
import math, os, sys
for d in reversed(os.environ.get("CRONOS_PYPATH", "").split(":")):
    if d and d not in sys.path:
        sys.path.insert(0, d)
from pymcpp import FFGraph, FFPartial, FFIntegral, exp as ffexp, sqrt as ffsqrt
from cronos import FFDom, FFModel, DAESLV, OCFESLV, CollocationOptions

npass = nfail = 0
def check(cond, what):
    global npass, nfail
    print("  %s  %s" % ("PASS" if cond else "FAIL", what))
    if cond: npass += 1
    else: nfail += 1
def rel(a, b):
    return abs(a - b) / max(1., abs(b))

OpP, OpI = FFPartial(), FFIntegral()
ALPHA, TEND = 0.1, 1.
TI = FFDom.ALL - FFDom.LB
ZI = FFDom.ALL - FFDom.LB - FFDom.UB
Role = FFModel.EqnRole
io, ib, ii = FFModel.EqnOptions(Role.INTERIOR), FFModel.EqnOptions(Role.BOUNDARY), FFModel.EqnOptions(Role.INITIAL)

def tight(S):
    o = S.options.SOLVE
    if isinstance(S, DAESLV):
        o.RTOL = o.ATOL = o.RTOLS = o.ATOLS = o.RTOLB = o.ATOLB = 1e-10
        o.NMAX = 20000
        o.VERBOSE = 0
    else:
        o.RES_TOL = 1e-12
    S.options.DISPLAY_LEVEL = 0

def rod(cls, inlet="param", nodes_t=3):
    """The rod u_t = ALPHA u_zz on z in [0,1], u(t,0) = inlet (a parameter p or a time-varying g(t)), u_z(t,1) = 0,
    u(0,z) = p (1-z)^2 (consistent with both boundary conditions, so the solution is smooth and J stays linear in
    p; with g(t) the initial field is 0); output J = int_0^1 u dz at TEND.  Returns (S, p_or_g, (t, z, u))."""
    G = FFGraph()
    t, z, u = G.add_var("t"), G.add_var("z"), G.add_var("u")
    S = cls(G)
    S.add_domain(t, FFDom(0., TEND, 2, FFDom.LGR, nodes_t))
    S.add_domain(z, FFDom(0., 1., 3, FFDom.LGL, 5))
    S.set_evolution_domain(t)
    S.add_state(u, [t, z], ref=0.5)
    if inlet == "param":
        w = G.add_var("p"); S.add_input(w, ref=1.)
    else:
        w = G.add_var("g"); S.add_input(w, [t], ref=1., type=FFDom.LGL, n_node=3)   # 3 levels per stage
    S.add_equation(OpP(u, t) - OpP(ALPHA * OpP(u, z), z), [t, z], [TI, ZI], io)
    S.add_equation(u - w,     [t, z], [TI, FFDom.LB], ib)
    S.add_equation(OpP(u, z), [t, z], [TI, FFDom.UB], ib)
    S.add_equation(u - (w * (1. - z) * (1. - z) if inlet == "param" else 0.), [t, z], [FFDom.LB, FFDom.ALL], ii)
    S.add_output(OpI(u, z), [t], point=[TEND])
    if cls is DAESLV:
        S.options.INTERFACE.IMPOSITION = CollocationOptions.IC_STRONG
        S.options.INTERFACE.SAT_SIGMA0 = 10.
    tight(S)
    return S, w, (t, z, u)

# ---- (A) Dirichlet inlet p: DAESLV == OCFESLV, forward == adjoint == central difference
print("(A) a rod with a Dirichlet inlet p, on DAESLV and on OCFESLV")
P = 0.7
S, p, (t, z, u) = rod(DAESLV)
check(S.setup(), "DAESLV setup of the PDE  (%s)" % (S.extract_error() or "ok"))
check(S.solve({p: P}) == DAESLV.NORMAL, "solve({p: P})")
J = S.val_function()[0]
O, po, _ = rod(OCFESLV, nodes_t=6)                      # OCFESLV at time order 6: the gap is its time error
check(O.setup(), "OCFESLV setup of the same declaration")
check(O.solve({po: P}) == OCFESLV.NORMAL, "OCFESLV solve({p: P})")
Jo = O.val_function()[0]
print("  J: DAESLV %.6e  OCFESLV %.6e  rel %.2e" % (J, Jo, rel(J, Jo)))
check(rel(J, Jo) < 1e-4, "DAESLV's output = OCFESLV's at time order 6, to 1e-4")
check(S.solve_fsens({p: P}) == DAESLV.NORMAL, "solve_fsens({p: P})")
gf = S.val_function_gradient()[0][0]                     # one row per parameter, one column per output
check(S.solve_asens({p: P}) == DAESLV.NORMAL, "solve_asens({p: P})")
ga = S.val_function_gradient()[0][0]
h = 1e-4
S.solve({p: P + h}); Jp = S.val_function()[0]
S.solve({p: P - h}); Jm = S.val_function()[0]
fd = (Jp - Jm) / (2 * h)
print("  dJ/dp: fwd %.6e  adj %.6e  FD %.6e" % (gf, ga, fd))
check(rel(gf, ga) < 1e-8, "forward = adjoint, to 1e-8")
check(rel(gf, fd) < 1e-6, "forward = central difference, to 1e-6")
# J is linear in p: dJ/dp = J/p exactly
check(rel(gf, J / P) < 1e-8, "dJ/dp = J/p (J linear in p), to 1e-8")

# ---- (B) a time-varying inlet g(t): one gradient row per DOF, each = a central difference
print("(B) a time-varying Dirichlet inlet g(t), 3 LGL levels on each of 2 stages")
S, g, (t, z, u) = rod(DAESLV, inlet="tv")
check(S.setup(), "setup  (%s)" % (S.extract_error() or "ok"))
dofs = S.input_dofs(g)
check(len(dofs) == 6, "6 declared DOFs for g (3 levels x 2 stages): %d" % len(dofs))
vals = [1. + 0.3 * k for k in range(len(dofs))]
check(S.solve_fsens({g: vals}) == DAESLV.NORMAL, "solve_fsens({g: [values by DOF]})")
J0 = S.val_function()[0]
gf = [r[0] for r in S.val_function_gradient()]            # one row per declared DOF of g, in input_dofs order
check(S.solve_asens({g: vals}) == DAESLV.NORMAL, "solve_asens({g: [values by DOF]})")
ga = [r[0] for r in S.val_function_gradient()]
check(len(gf) == len(dofs) and len(ga) == len(dofs), "one gradient row per declared DOF, forward and adjoint")
fd, h = [], 1e-4
for k in range(len(dofs)):
    vp = list(vals); vp[k] += h; S.solve({g: vp}); Jp = S.val_function()[0]
    vm = list(vals); vm[k] -= h; S.solve({g: vm}); Jm = S.val_function()[0]
    fd.append((Jp - Jm) / (2 * h))
scale = max(1e-12, max(abs(v) for v in fd))
e_fa = max(abs(a - b) for a, b in zip(gf, ga)) / scale
e_fd = max(abs(a - b) for a, b in zip(gf, fd)) / scale
for k in range(len(dofs)):
    print("  DOF %d: fwd %.4e  adj %.4e  FD %.4e" % (k, gf[k], ga[k], fd[k]))
print("  max scaled: adj-fwd %.2e  fwd-FD %.2e" % (e_fa, e_fd))
check(e_fa < 1e-8, "adjoint = forward, to 1e-8 scaled by the largest gradient")
check(e_fd < 1e-5, "forward = central differences, to 1e-5 scaled")

# ---- (C) a transition on the field: closed forms exact on any mesh
print("(C) a transition u(tau+) = u(tau-)/2 at tau = 0.5 on a zero-flux rod with a source p e^-t")
TAU, PS = 0.5, 0.7
eT, eTau = math.exp(-TEND), math.exp(-TAU)
Jex  = 0.5 * (1. + PS * (1. - eTau)) + PS * (eTau - eT)
dJex = 0.5 * (1. - eTau) + eTau - eT
def rod_tr():
    G = FFGraph()
    t, z, u, p = G.add_var("t"), G.add_var("z"), G.add_var("u"), G.add_var("p")
    S = DAESLV(G)
    S.add_domain(t, FFDom(0., TEND, 2, FFDom.LGR, 3))
    S.add_domain(z, FFDom(0., 1., 3, FFDom.LGL, 5))
    S.set_evolution_domain(t)
    S.add_state(u, [t, z], ref=1.)
    S.add_input(p, ref=1.)
    S.add_equation(OpP(u, t) - OpP(ALPHA * OpP(u, z), z) - p * ffexp(-t), [t, z], [TI, ZI], io)
    S.add_equation(OpP(u, z), [t, z], [TI, FFDom.LB], ib)
    S.add_equation(OpP(u, z), [t, z], [TI, FFDom.UB], ib)
    S.add_equation(u - 1.,    [t, z], [FFDom.LB, FFDom.ALL], ii)
    S.add_transition([0.5 * u], [u], t, TAU)
    S.add_output(OpI(u, z), [t], point=[TEND])
    S.options.INTERFACE.IMPOSITION = CollocationOptions.IC_STRONG
    S.options.INTERFACE.SAT_SIGMA0 = 10.
    tight(S)
    return S, p
S, p = rod_tr()
check(S.setup(), "setup with a transition  (%s)" % (S.extract_error() or "ok"))
check(S.solve({p: PS}) == DAESLV.NORMAL, "solve")
J = S.val_function()[0]
check(S.solve_fsens({p: PS}) == DAESLV.NORMAL, "solve_fsens across the stop")
gf = S.val_function_gradient()[0][0]
check(S.solve_asens({p: PS}) == DAESLV.NORMAL, "solve_asens across the stop")
ga = S.val_function_gradient()[0][0]
print("  J %.8e (closed form %.8e)   dJ/dp fwd %.8e adj %.8e (closed form %.8e)" % (J, Jex, gf, ga, dJex))
check(rel(J, Jex) < 1e-8, "J = closed form, to 1e-8")
check(rel(gf, dJex) < 1e-8 and rel(ga, dJex) < 1e-8, "dJ/dp = closed form, forward and adjoint, to 1e-8")

# ---- (D) the semi-discrete system
print("(D) the semi-discrete system's unknowns and the field's nodes")
S, p, (t, z, u) = rod(DAESLV)
check(S.setup(), "setup  (%s)" % (S.extract_error() or "ok"))
xc = S.semidiscrete_states()
nodes = S.semidiscrete_nodes(u)
print("  %d semi-discrete unknowns; u has %d nodes" % (len(xc), len(nodes)))
check(len(xc) == S.nx() and len(xc) > 0, "semidiscrete_states has nx entries")
check(len(nodes) == 15, "u has 15 nodes on 3 LGL(5) elements (each element owns its end nodes): %d" % len(nodes))
cols = [c for c, _ in nodes]
zs = [coord[z] for _, coord in nodes]
check(all(0 <= c < len(xc) for c in cols) and len(set(cols)) == len(cols), "node columns valid and distinct")
check(all(zs[i+1] >= zs[i] - 1e-12 for i in range(len(zs) - 1)) and abs(zs[0]) < 1e-14 and abs(zs[-1] - 1.) < 1e-14,
      "node coordinates non-decreasing from 0 to 1 (mesh order, the interface nodes doubled)")

# ---- (E) nonlinear functions of horizon integrals: on the zero-flux rod with u0 = 1 and source p e^-t the
#      inventory is J(t) = 1 + p (1 - e^-t) exactly, so A = int J dt, B = int J^2 dt, sqrt A and A/B are closed
print("(E) nonlinear functions of horizon integrals: A, sqrt A, A/B with B = int (int u dz)^2 dt")
PE = 0.7; e1, e2 = math.exp(-TEND), math.exp(-2. * TEND)
Ae, dAe = 1. + PE * e1, e1
Be = (1 + PE) ** 2 - 2 * PE * (1 + PE) * (1 - e1) + PE * PE * (1 - e2) / 2.
dBe = 2 * (1 + PE) - 2 * (1 + 2 * PE) * (1 - e1) + PE * (1 - e2)
Fe = [Ae, math.sqrt(Ae), Ae / Be]
dFe = [dAe, dAe / (2 * math.sqrt(Ae)), (dAe * Be - Ae * dBe) / (Be * Be)]
def rod_nl():
    G = FFGraph()
    t, z, u, p = G.add_var("t"), G.add_var("z"), G.add_var("u"), G.add_var("p")
    S = DAESLV(G)
    S.add_domain(t, FFDom(0., TEND, 2, FFDom.LGR, 3))
    S.add_domain(z, FFDom(0., 1., 3, FFDom.LGL, 5))
    S.set_evolution_domain(t)
    S.add_state(u, [t, z], ref=1.)
    S.add_input(p, ref=1.)
    S.add_equation(OpP(u, t) - OpP(ALPHA * OpP(u, z), z) - p * ffexp(-t), [t, z], [TI, ZI], io)
    S.add_equation(OpP(u, z), [t, z], [TI, FFDom.LB], ib)
    S.add_equation(OpP(u, z), [t, z], [TI, FFDom.UB], ib)
    S.add_equation(u - 1.,    [t, z], [FFDom.LB, FFDom.ALL], ii)
    J = OpI(u, z); A = OpI(J, t); B = OpI(J * J, t)
    S.add_output(A); S.add_output(ffsqrt(A)); S.add_output(A / B)
    S.options.INTERFACE.IMPOSITION = CollocationOptions.IC_STRONG
    S.options.INTERFACE.SAT_SIGMA0 = 10.
    tight(S)
    return S, p
S, p = rod_nl()
check(S.setup() and S.nf() == 3, "three nonlinear outputs set up (nf = 3)  (%s)" % (S.extract_error() or "ok"))
check(S.solve({p: PE}) == DAESLV.NORMAL, "solve")
F = list(S.val_function())
check(S.solve_fsens({p: PE}) == DAESLV.NORMAL, "solve_fsens")
dF = list(S.val_function_gradient()[0])
check(S.solve_asens({p: PE}) == DAESLV.NORMAL, "solve_asens")
dA = list(S.val_function_gradient()[0])
print("  A %.8e sqrtA %.8e A/B %.8e   d/dp fwd %.6e %.6e %.6e" % (F[0], F[1], F[2], dF[0], dF[1], dF[2]))
check(max(rel(a, b) for a, b in zip(F, Fe)) < 1e-8, "values = closed forms, to 1e-8")
check(max(rel(a, b) for a, b in zip(dF, dFe)) < 1e-8, "forward gradients = closed forms, to 1e-8")
check(max(rel(a, b) for a, b in zip(dA, dFe)) < 1e-6, "adjoint gradients = closed forms, to 1e-6")

print("\ncheck_daeslv_mol: %d passed, %d failed -- %s" % (npass, nfail, "FAILURES" if nfail else "ALL PASS"))
sys.exit(1 if nfail else 0)
