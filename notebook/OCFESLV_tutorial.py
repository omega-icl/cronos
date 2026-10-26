#!/usr/bin/env python
# coding: utf-8

# # Nonlinear Heat Conduction in a Rod with a Heater Schedule
# 
# This tutorial solves a nonlinear parabolic PDE with the orthogonal-collocation-on-finite-elements solver `OCFESLV` of the library `cronos` (CRONOS), the model being declared on a DAG of `pymcpp` (MC++). The PDE is solved both **monolithically** (all time elements at once) and by **marching** (one time element after the other); the solution field is read back anywhere through `eval_solution`; and the outputs are differentiated with respect to a time-varying control.
# 
# A rod $z \in [0,1]$ is held at temperature $0$ at $z=0$ and heated at $z=1$ by a heater of power $q(t)$. Its thermal conductivity increases with temperature, so heat conduction is **nonlinear**:
# $$\begin{align*}
#   \frac{\partial u}{\partial t} &= \frac{\partial}{\partial z}\left( D(u)\,\frac{\partial u}{\partial z} \right), \qquad D(u) = D_0\,(1 + \beta u), && (t,z) \in (0,2]\times (0,1) \\
#   u(t,0) &= 0, \qquad D(u)\,\frac{\partial u}{\partial z}(t,1) = q(t) && t \in (0,2]\\
#   u(0,z) &= a\,z && z \in [0,1]
# \end{align*}$$
# with $D_0 = 0.1$ and $\beta = 1$. 
# 
# The initial profile is the linear one that carries the initial heater power, $D_0(1+\beta a)\,a = q(0)$, so that the initial and boundary data are compatible.
# 
# The heater is cycled, from full power down to a quarter of it and back, twice over the horizon: $q(t) = 0.625 + 0.375\cos(2\pi t)$.
# 
# The outputs of interest are the temperature at the heated end at the final time, $u(2,1)$, and the space-time mean temperature $\int_0^2\!\int_0^1 u(t,z)\,dt\,dz$.

# We start by importing the `pymcpp` and `cronos` libraries.

# In[1]:


import pymcpp
import cronos
from pymcpp import FFGraph, FFPartial, FFEval, FFIntegral
from cronos import FFDom, FFModel, OCFESLV, FFOCFESLV


# ## The model
# 
# The model lives on a DAG: the independent variables $t$ and $z$, the state $u(t,z)$ and the input $q(t)$. Derivatives, point evaluations and integrals along the domains are DAG operations: `FFPartial`, `FFEval` and `FFIntegral`.
# 
# Each domain is partitioned into finite elements, with collocation nodes of a given family in each element: Legendre-Gauss-Radau (`LGR`) in time, which suits an evolution direction, and Legendre-Gauss-Lobatto (`LGL`) in space, whose nodes include the element ends.

# In[2]:


import numpy as np

# initial profile a z, with D0 (1 + beta a) a = q(0) = 1
D0, beta = 0.1, 1.0
a = (-1. + np.sqrt(1. + 4. * beta / D0)) / (2. * beta)
print("initial slope a = %.6f" % a)

def heater(t):
    """Heater power at time t"""
    return 0.625 + 0.375 * np.cos(2. * np.pi * t)


# The solver `OCFESLV` is populated through the `FFModel` interface: the domains (with the evolution direction $t$), the state and the input with their reference values, the equations with their region and role, and the outputs.
# 
# The flux form $\partial_z(D(u)\,\partial_z u)$ is written as it stands, as the $z$-derivative of the DAG expression $D(u)\,\partial_z u$; the solver expands it and reduces the second derivative to first order automatically (an auxiliary state $\partial_z u$ is introduced at `setup()`). The regions are given by masks: `FFDom.ALL - FFDom.LB` is $(0,1]$, and so on. The heater input is discretised on the time elements with 3 Lobatto nodes each (`LGL`, `n_node=3`): a quadratic in every element, with nodes at both element ends, so that the heater schedule is continuous from one element to the next.

# In[3]:


G = FFGraph()
t, z = G.add_var("t"), G.add_var("z")
u, q = G.add_var("u"), G.add_var("q")

S = OCFESLV(G)
OpP, OpE, OpI = FFPartial(), FFEval(), FFIntegral()

NT, NZ = 10, 8                       # time and space elements
S.add_domain(t, FFDom(0., 2., NT, FFDom.LGR, 4))
S.add_domain(z, FFDom(0., 1., NZ, FFDom.LGL, 5))
S.set_evolution_domain(t)

S.add_state(u, [t, z], ref=0.)
S.add_input(q, [t], ref=1., type=FFDom.LGL, n_node=3)    # 3 nodes per time element, both ends included

D = D0 * (1. + beta * u)
S.add_equation(OpP(u, t) - OpP( D * OpP(u, z), z), [t, z], [FFDom.ALL - FFDom.LB, FFDom.ALL - FFDom.LB - FFDom.UB], FFModel.EqnRole.INTERIOR)
S.add_equation(u,                                  [t, z], [FFDom.ALL - FFDom.LB, FFDom.LB],                        FFModel.EqnRole.BOUNDARY)          # u(t,0) = 0
S.add_equation(D * OpP(u, z) - q,                  [t, z], [FFDom.ALL - FFDom.LB, FFDom.UB],                        FFModel.EqnRole.BOUNDARY)          # heater flux at z = 1
S.add_equation(u - a * z,                          [t, z], [FFDom.LB, FFDom.ALL],                                   FFModel.EqnRole.INITIAL)           # u(0,z) = a z

S.add_output(OpE(u, {t: 1, z: 1}, {t: 2., z: 1.}))        # u(2,1)
S.add_output(OpI(u, {t: 1, z: 1}))                        # space-time mean temperature

S.options.INTERFACE.IMPOSITION = OCFESLV.Options.IC_STRONG
S.options.SOLVE.MARCHING       = False
S.options.SOLVE.RES_TOL        = 1e-10


# Interface continuity between elements is imposed **strongly** (`IC_STRONG`): exact continuity rows. The alternative, weak imposition through penalties (`IC_WEAK`, the current default), depends on a penalty strength and is markedly less accurate on this problem at this resolution.

# ## Monolithic solution
# 
# `setup()` analyses the model and builds the collocation system; the model report shows what the solver sees: the equations written out with their regions, the auxiliary state introduced by the order reduction, the classification and the degrees of freedom.

# In[4]:


S.options.DISPLAY_LEVEL = 1
S.setup()


# The unknowns live in NumPy arrays: `init()` returns the collocation state vector `var`, initialised from the reference values, and the input vector `inp`. `sample_input` writes the heater schedule into `inp` -- a function of the coordinates, evaluated at the input's collocation nodes -- and `solve` overwrites `var` with the solution:

# In[5]:


import time

var, inp = S.init()
S.sample_input(q, lambda c: heater(c[t]), inp)
print( inp )

t0 = time.perf_counter()
report = S.solve(var, inp)
cpu_mono = time.perf_counter() - t0
f_mono = S.val_functions()
print(report, " (%.2f s)" % cpu_mono)
print("u(2,1) = %.8f   u average = %.8f" % tuple(f_mono))


# ## The solution field
# 
# `eval_solution(state, point)` reads the converged solution at any point of the domain -- here on a regular grid, for a map of the temperature field:

# In[6]:


import matplotlib.pyplot as plt

tt = np.linspace(0., 2., 161)
zz = np.linspace(0., 1., 41)
U = np.array([[S.eval_solution(u, {t: tv, z: zv}) for tv in tt] for zv in zz])

fig, ax = plt.subplots(figsize=(9, 4))
pc = ax.pcolormesh(tt, zz, U, shading="auto", cmap="viridis")
ax.set_xlabel("$t$")
ax.set_ylabel("$z$")
fig.colorbar(pc, label="$u(t,z)$")
plt.show()


# The same field as a surface over $(t, z)$: the rod warms from the heated end while the heater is near full power and cools as it is turned down, twice over the horizon, and the temperature at $z = 0$ stays pinned at zero.
# 

# In[7]:


from matplotlib import cm

TT, ZZ = np.meshgrid(tt, zz)
fig = plt.figure(figsize=(11, 6))
ax = fig.add_subplot(projection="3d")
surf = ax.plot_surface(TT, ZZ, U, cmap=cm.viridis, rstride=1, cstride=2, linewidth=0.2, edgecolor="k", alpha=0.95)
ax.set_xlabel("$t$")
ax.set_ylabel("$z$")
ax.set_zlabel("$u(t,z)$")
ax.view_init(elev=28, azim=-130)
fig.colorbar(surf, shrink=0.6, pad=0.1, label="$u(t,z)$")
fig.subplots_adjust(left=0.08)
plt.show()


# ## Marching
# 
# For an evolution problem, the collocation system can also be solved one time element after the other (`options.SOLVE.MARCHING`), each window starting from the end of the previous one. This is much cheaper: each solve is small. With `options.OUTPUT.MARCH_STORE` (on by default), every window is kept, so `eval_solution` reads the full marched field as well:

# In[8]:


S.options.SOLVE.MARCHING  = True
S.options.DISPLAY_LEVEL   = 1
S.setup()

var, inp = S.init()
S.sample_input(q, lambda c: heater(c[t]), inp)
#print( inp )

t0 = time.perf_counter()
report = S.solve(var, inp)
cpu_march = time.perf_counter() - t0
f_march = S.val_functions()
print(report, " (%.2f s, %d windows)" % (cpu_march, S.n_march_steps()))
print("u(2,1) = %.8f   u average = %.8f" % tuple(f_march))

U2 = np.array([[S.eval_solution(u, {t: tv, z: zv}) for tv in tt] for zv in zz])
print("monolithic: %.3f s   marching: %.3f s   max |u_mono - u_march| over the grid: %.2e"
      % (cpu_mono, cpu_march, np.abs(U - U2).max()))


# Marching is causal -- each window depends only on the past -- and so is the monolithic time discretisation: every time element starts from the previous element's end state, as a marching window does. The two modes therefore solve the same discretisation and agree to round-off. Marching is several times faster.

# ## Sensitivities with respect to the heater schedule
# 
# The heater input is registered as a **control**: its degrees of freedom are its values at the collocation nodes, 3 per time element. `solve_fsens` (forward) and `solve_asens` (adjoint) then compute the derivatives of the outputs with respect to them, returned by `sens_jacobian()` as an array of shape (outputs, control DOFs):

# In[9]:


S.options.DISPLAY_LEVEL = 0
S.register_control(q)
S.setup()

var, inp = S.init()
S.sample_input(q, lambda c: heater(c[t]), inp)

S.solve_fsens(var.copy(), inp)
J_fwd = S.sens_jacobian()
S.solve_asens(var.copy(), inp)
J_adj = S.sens_jacobian()
print("sens_jacobian:", J_fwd.shape, "| max |forward - adjoint| = %.2e" % np.abs(J_fwd - J_adj).max())


# A finite-difference check on one degree of freedom: `encode_controls` reads the control values out of `inp`, `decode_controls` writes them back:

# In[10]:


S.solve(var.copy(), inp)
f0 = np.array(S.val_functions())

k, h = -1, 1e-6
p = S.encode_controls(inp)
p[k] += h
inp_h = inp.copy()
S.decode_controls(p, inp_h)

S.solve(var.copy(), inp_h)
fd = (np.array(S.val_functions()) - f0) / h
print("d outputs / d q_%d: forward %s   finite difference %s" % (k, np.round(J_fwd[:, k], 6), np.round(fd, 6)))


# ## The heated rod as a DAG
# 
# The solver can also be embedded in a DAG as an external operation, with `FFOCFESLV`: its outputs become DAG variables, functions of the DAG variables given for its inputs, which can then be composed with other DAG operations, evaluated, sampled and differentiated like any other DAG function.
# 
# Here the heater schedule is parametrised in the DAG as $q(t) = A + B\cos(2\pi t)$, with $A$ and $B$ DAG variables; $A = 0.625$, $B = 0.375$ is the schedule used so far. The heater's degrees of freedom -- its values at the time nodes `t_dof` found above -- become DAG expressions $A + B\cos(2\pi t_k)$, passed in the order of the control DOFs. We embed a **marching** solver, the faster of the two. Under the default policy `SHALLOW` the operation refers to the solver, which the DAG keeps alive; `COPY` gives the operation its own copy -- the policy to use when embedding the same solver several times with different inputs.

# In[17]:


DAG = FFGraph()
A = DAG.add_var("A")
B = DAG.add_var("B")

S.options.MAXTHREAD = 1        # one BLAS/OpenMP thread per solve: veval below runs the solves themselves in parallel
S.setup()

t_grid = inp.copy()
S.sample_input(q, lambda c: c[t], t_grid)               # the input's own coordinates, as an input vector
t_dof = np.array(S.encode_controls(t_grid))             # the time of each heater DOF
Q = [A + B * np.cos(2. * np.pi * tk) for tk in t_dof]   # the heater DOFs as DAG expressions, at the times found above

OpOCFE = FFOCFESLV()
[u21, ubar] = OpOCFE({q: Q}, S)

u21.set("u21")
ubar.set("ubar")
#DAG.output([u21, ubar])
print(u21, " = ", u21.str())


# A derived output, the ratio of the final temperature at the heated end to the space-time mean:

# In[12]:


uratio = u21 / ubar
uratio.set("uratio")


# Evaluated at the nominal schedule, the DAG agrees with the solves above:

# In[13]:


[u21_val, ubar_val, uratio_val] = DAG.eval([u21, ubar, uratio], [A, B], [0.625, 0.375])
print("u(2,1) = %.8f   mean = %.8f   ratio = %.6f" % (u21_val, ubar_val, uratio_val))
print("from the solve above: u(2,1) = %.8f   mean = %.8f" % tuple(f0))


# Multiple DAG evaluations over the schedule parameters $(A, B) \in [0.25, 1] \times [0, 0.5]$; each evaluation is a complete simulation of the rod:
# 
# With `DAG.options.MAXTHREAD = 0`, `veval` spreads the evaluations over all available threads, each thread with its own copy of the embedded solver. The solver itself was therefore limited to a single BLAS/OpenMP thread above (`S.options.MAXTHREAD = 1`): otherwise each concurrent solve would also start its own pool of linear-algebra threads, oversubscribing the cores.

# In[14]:


from scipy.stats import qmc

sampler = qmc.Sobol(d=2, scramble=False, seed=1)
samInput = sampler.random_base2(m=7)                  # 2^7 = 128 evaluations
samInput[:, 0] = 0.25 + samInput[:, 0] * (1. - 0.25)
samInput[:, 1] = samInput[:, 1] * 0.5

DAG.options.MAXTHREAD = 0                             # use all available threads
samOutput = np.array(DAG.veval([u21, uratio], [A, B], samInput.tolist(), walltime=True))
print(samOutput[:5])


# In[15]:


fig = plt.figure(figsize=(12, 5))
ax = fig.add_subplot(1, 2, 1, projection='3d')
ax.scatter(samInput[:, 0], samInput[:, 1], samOutput[:, 0], s=4)
ax.set_xlabel("$A$"); ax.set_ylabel("$B$"); ax.set_zlabel("$u(2,1)$")
ax.view_init(elev=25, azim=230)
ax2 = fig.add_subplot(1, 2, 2)
sc = ax2.scatter(samInput[:, 0], samInput[:, 1], c=samOutput[:, 1], s=18, cmap="viridis")
ax2.set_xlabel("$A$"); ax2.set_ylabel("$B$"); ax2.set_title(r"$u(2,1)\,/\,\bar u$")
fig.colorbar(sc, ax=ax2)
plt.tight_layout()
plt.show()


# The derivatives of the outputs with respect to $A$ and $B$ come from the DAG (`fdiff`), which differentiates the embedded solve through the solver's sensitivities. They agree with the chain rule applied to the sensitivity Jacobian computed earlier, since $\partial q_k/\partial A = 1$ and $\partial q_k/\partial B = \cos(2\pi t_k)$:

# In[16]:


[Row, Col, Grad] = DAG.fdiff([u21, ubar], [A, B])
Grad_val = DAG.eval(Grad, [A, B], [0.625, 0.375])
J_AB = np.zeros((2, 2))
for i, j, g in zip(Row, Col, Grad_val):
    J_AB[i, j] = g
print("DAG fdiff:\n", J_AB)

J_chain = np.column_stack([J_fwd.sum(axis=1), J_fwd @ np.cos(2. * np.pi * t_dof)])
print("chain rule on the sensitivity Jacobian:\n", J_chain)
print("max |DAG fdiff - chain rule| = %.1e" % np.abs(J_AB - J_chain).max())

