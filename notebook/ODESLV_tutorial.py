# %% [markdown]
# # Simulation of a Serial Batch Reactor Network
#
# The purpose of this tutorial is to implement and simulate a serial batch chemical reactor network using the libraries `pymcpp` (MC++) and `cronos` (CRONOS).
#
# The network consists of two batch reactors connected in series. The reaction mechanism characterizing the two reactors is the same:
# $$\begin{equation*}
#     2{\sf A} \xrightarrow{k_1} {\sf B} \xrightarrow{k_2} {\sf C}
# \end{equation*}$$
# where component ${\sf B}$ is the desired product.
#
# The reactor dynamics are ordinary differential equations (ODEs) that describe the evolution of component molar concentrations in continuous time:
# $$\begin{align*}
#         \dot{c}_{{\sf A},j}(t) &= -2k_{1,j}(T_j)\,c_{{\sf A},j}(t)^2\\
#         \dot{c}_{{\sf B},j}(t) &= k_{1,j}(T_j)\,c_{{\sf A},j}(t)^2 - k_{2,j}(T_j)\,c_{{\sf B},j}(t)\\
#         \dot{c}_{{\sf C},j}(t) &= k_{2,j}(T_j)\,c_{{\sf B},j}(t)
# \end{align*}$$
# where $c_{i,j}(t)$ ($\text{kmol m}^{-3}$) describes the molar concentration of component $i={\sf A},{\sf B},{\sf C}$ in reactor $j=1,2$; and $T_j$ (K) the temperature in reactor $j$.
#
# The system adheres to Arrhenius kinetics:
# $$\begin{equation*}
#     k_{r,j}(T) = k^\circ_{r,j} \exp\left(-\frac{E_{r,j}}{RT}\right)
# \end{equation*}$$
# where $T$ (K) is the reaction temperature; $k_{r,j}$ the temperature-dependent rate constant of reaction step $r=1,2$ for reactor $j$; $k^\circ_{r,j}$ the pre-exponential factor; $E_{r,j}$ the activation energy; and $R$ denotes the ideal gas constant.
#
# The time is normalized by the batch duration $\tau_j$, $t \in [0,1]$, so that $\tau_j$ appears as a parameter of the dynamics.

# %% [markdown]
# We start by importing both the `pymcpp` and `cronos` libraries. Set `CRONOS_PYPATH` (or edit the cell) if the modules are not on the Python path.

# %%
import pymcpp
import cronos
from pymcpp import FFGraph, FFPartial, FFEval
from cronos import FFDom, FFModel, ODESLV, FFODESLV

# %% [markdown]
# ## A single batch reactor
#
# The model is declared on a DAG: the independent variable $t$, the states, the parameters (batch duration, temperature and initial concentrations) and the constants (kinetic parameters). The derivatives and point evaluations along $t$ are DAG operations: `FFPartial` and `FFEval`.

# %%
IVPDAG = FFGraph()

# Independent variable (normalized time)
t = IVPDAG.add_var("t")

# States
cA = IVPDAG.add_var("cA")
cB = IVPDAG.add_var("cB")
cC = IVPDAG.add_var("cC")

# Parameters
tau = IVPDAG.add_var("tau")
T   = IVPDAG.add_var("T")
cA0 = IVPDAG.add_var("cA0")
cB0 = IVPDAG.add_var("cB0")
cC0 = IVPDAG.add_var("cC0")

# Constants
kr1 = IVPDAG.add_var("kr1")
Ea1 = IVPDAG.add_var("Ea1")
kr2 = IVPDAG.add_var("kr2")
Ea2 = IVPDAG.add_var("Ea2")

OpP, OpE = FFPartial(), FFEval()

# %% [markdown]
# A solver `ODESLV` is created and populated with the parametric initial value problem, through the `FFModel` interface: the time domain, the states, the inputs (parameters) and constants, the equations with their region and role, and the outputs:

# %%
IVP = ODESLV(IVPDAG)
IVP.add_domain(t, FFDom(0., 1., 1))    # normalized time horizon, one stage
IVP.set_evolution_domain(t)
IVP.add_state(cA, [t])
IVP.add_state(cB, [t])
IVP.add_state(cC, [t])
for p in (tau, T, cA0, cB0, cC0):
    IVP.add_input(p)
IVP.set_constant([kr1, Ea1, kr2, Ea2])

R = 8.314  # J/mol·K
k1 = kr1 * pymcpp.exp(-Ea1 / (R * T))
k2 = kr2 * pymcpp.exp(-Ea2 / (R * T))

INTERIOR = FFModel.EqnOptions(FFModel.EqnRole.INTERIOR)
INITIAL  = FFModel.EqnOptions(FFModel.EqnRole.INITIAL)
T_INT = FFDom.ALL - FFDom.LB           # (0,1]: the dynamics
IVP.add_equation(OpP(cA, t) - tau * (-2 * k1 * cA**2),         [t], [T_INT], INTERIOR)
IVP.add_equation(OpP(cB, t) - tau * (k1 * cA**2 - k2 * cB),    [t], [T_INT], INTERIOR)
IVP.add_equation(OpP(cC, t) - tau * (k2 * cB),                 [t], [T_INT], INTERIOR)
IVP.add_equation([cA - cA0, cB - cB0, cC - cC0], [t], [FFDom.LB], INITIAL)   # t = 0: initial values

for c in (cA, cB, cC):                 # batch-end concentrations
    IVP.add_output(OpE(c, t, 1.))

IVP.setup()
IVP.report()

# %% [markdown]
# The model report shows the declared model as the solver sees it: the equations written out with their regions, the outputs, the classification and the degrees of freedom.

# %% [markdown]
# A simulation of the first reactor for a batch duration $\tau=600~\rm min$ and temperature $T=750~\rm K$, starting from $c_{\sf A}=2~\text{kmol m}^{-3}$. The inputs and constants are given **by name**; the forward (`solve_fsens`) and adjoint (`solve_asens`) sensitivity analyses give the gradients of the batch-end concentrations with respect to the parameters:

# %%
values = {tau: [600.], T: [750.], cA0: [2.], cB0: [0.], cC0: [0.],
          kr1: [6.66e-3], Ea1: [2.52e3], kr2: [1.03], Ea2: [5e3]}

IVP.options.DISPLAY = 1      # displays numerical integration results
IVP.options.RESRECORD = 200  # record 200 points along time horizon

IVP.solve(values)
IVP.solve_fsens(values)
IVP.solve_asens(values)

# %% [markdown]
# The options of `ODESLV` combine the model options (`FFModel.Options`) and the CVODES options; each is documented (e.g. `help(ODESLV.Options.RTOL)`, or `help(IVP.options)` for all of them):

# %%
for name in ("INTMETH", "NLINSOL", "LINSOL", "RTOL", "ATOL", "NMAX", "RESRECORD", "DISPLAY"):
    print("%-10s = %-24s %s" % (name, getattr(IVP.options, name), getattr(ODESLV.Options, name).__doc__.split(".")[0]))

# %% [markdown]
# The gradients of the batch-end concentrations, one row per parameter and one column per output, agree between forward and adjoint sensitivity analysis:

# %%
import numpy as np
import matplotlib.pyplot as plt

IVP.options.DISPLAY = 0

IVP.solve_fsens(values)
G_fwd = np.array(IVP.val_function_gradient())

IVP.solve_asens(values)
G_adj = np.array(IVP.val_function_gradient())

names = [str(p) for p in IVP.var_parameter]
print("%-5s %14s %14s %14s" % ("", "d cA(1)", "d cB(1)", "d cC(1)"))
for n, row in zip(names, G_fwd):
    print("%-5s %14.6e %14.6e %14.6e" % ((n,) + tuple(row)))
print("max |forward - adjoint| =", np.abs(G_fwd - G_adj).max())

# %% [markdown]
# ### Trajectories
#
# The trajectories are displayed from the `results_solve` field: NumPy arrays of the recorded time points `t`, states `x` (in the order of `var_state`) and quadratures `q`. The normalized time is scaled back by $\tau$:

# %%
IVP.solve(values)

TAU = values[tau][0]
res = IVP.results_solve
print("t:", res.t.shape, "| x:", res.x.shape)
#print( res.x )

plt.xlabel('$t$ (min)')
plt.ylabel('$c$ (kmol m$^{-3}$)')
plt.plot(res.t * TAU, res.x[:, 0], "-b", label="$c_{A,1}$")
plt.plot(res.t * TAU, res.x[:, 1], "-r", label="$c_{B,1}$")
plt.plot(res.t * TAU, res.x[:, 2], "-g", label="$c_{C,1}$")
plt.legend(loc="best")
plt.show()

# %% [markdown]
# The sensitivity trajectories of `solve_fsens` are in `results_fsens`: `xp[i]` holds the sensitivities of the states with respect to the $i$-th sensitivity direction (here the parameters, in the order of `var_parameter`). The sensitivities of the concentrations with respect to the temperature $T$:

# %%
IVP.solve_fsens(values)
sens = IVP.results_fsens
iT = names.index("T")
print("xp:", sens.xp.shape, "(directions, points, states)")

plt.xlabel('$t$ (min)')
plt.ylabel(r'$\partial c / \partial T$ (kmol m$^{-3}$ K$^{-1}$)')
plt.plot(sens.t * TAU, sens.xp[iT, :, 0], "-b", label=r"$\partial c_{A,1}/\partial T$")
plt.plot(sens.t * TAU, sens.xp[iT, :, 1], "-r", label=r"$\partial c_{B,1}/\partial T$")
plt.plot(sens.t * TAU, sens.xp[iT, :, 2], "-g", label=r"$\partial c_{C,1}/\partial T$")
plt.legend(loc="best")
plt.show()

# %%
IVP.solve_asens(values)
adj = IVP.results_asens
iF = 2 # Output cC(1)
print("l:", adj.l.shape, "(directions, points, states)")

plt.xlabel('$t$ (min)')
plt.ylabel(r'$\lambda$')
plt.plot(adj.t * TAU, adj.l[iT, :, 0], "-b", label=r"$\lambda_{A,1}$")
plt.plot(adj.t * TAU, adj.l[iT, :, 1], "-r", label=r"$\lambda_{B,1}$")
plt.plot(adj.t * TAU, adj.l[iT, :, 2], "-g", label=r"$\lambda_{C,1}$")
plt.legend(loc="best")
plt.show()

# %% [markdown]
# ### Both reactors in series
#
# The second reactor starts from the batch-end concentrations of ${\sf A}$ and ${\sf B}$ in the first one, without ${\sf C}$ (separated between the reactors). Two successive simulations with the same solver give the trajectories through the network, here for $\tau_1=\tau_2=400~\rm min$ and $T_1=T_2=800~\rm K$:

# %%
kin = {kr1: [6.66e-3], Ea1: [2.52e3], kr2: [1.03], Ea2: [5e3]}
TAU1, T1v, TAU2, T2v = 400., 800., 400., 800.

IVP.solve({tau: [TAU1], T: [T1v], cA0: [2.], cB0: [0.], cC0: [0.], **kin})
r1 = IVP.results_solve                  # NumPy arrays: copies, kept after the next solve
cA1f, cB1f, cC1f = IVP.val_function()

IVP.solve({tau: [TAU2], T: [T2v], cA0: [cA1f], cB0: [cB1f], cC0: [0.], **kin})
r2 = IVP.results_solve
traj1 = np.column_stack([r1.t * TAU1, r1.x])
traj2 = np.column_stack([TAU1 + r2.t * TAU2, r2.x])

plt.xlabel('$t$ (min)')
plt.ylabel('$c$ (kmol m$^{-3}$)')
for k, (col, lab) in enumerate((("b", "A"), ("r", "B"), ("g", "C")), start=1):
    plt.plot(traj1[:, 0], traj1[:, k], "-" + col, label="$c_{%s}$" % lab)
    plt.plot(traj2[:, 0], traj2[:, k], "--" + col)
plt.axvline(TAU1, color="grey", lw=0.8)
plt.text(TAU1 / 2, plt.ylim()[1] * 0.95, "reactor 1", ha="center")
plt.text(TAU1 + TAU2 / 2, plt.ylim()[1] * 0.95, "reactor 2", ha="center")
plt.legend(loc="center right")
plt.show()

# %% [markdown]
# ## The reactor network as a DAG
#
# To simulate the two batch reactors connected in series as one function, we define a separate DAG and embed the solver twice as an external operation `FFODESLV`. Each embedding maps the solver's inputs to DAG variables: map 1 lists the inputs differentiated through the solver's sensitivities, map 2 the remaining ones (here the constants), whose values pass through the DAG:

# %%
DAG = FFGraph()
OpIVP = FFODESLV()
IVP.options.DISPLAY = 0    # turn off display during numerical simulation
IVP.options.RESRECORD = 0  # turn off trajectory record

# %%
# reactor #1 initial concentrations
cA1_0 = DAG.add_var("cA1_0")
cB1_0 = DAG.add_var("cB1_0")
cC1_0 = DAG.add_var("cC1_0")

# reactor #1 parameters
tau1 = DAG.add_var("tau1")
T1   = DAG.add_var("T1")

# reactor #1 constants
kr11 = DAG.add_var("kr11")
Ea11 = DAG.add_var("Ea11")
kr21 = DAG.add_var("kr21")
Ea21 = DAG.add_var("Ea21")

[cA1_f, cB1_f, cC1_f] = OpIVP({tau: [tau1], T: [T1], cA0: [cA1_0], cB0: [cB1_0], cC0: [cC1_0]},
                              {kr1: [kr11], Ea1: [Ea11], kr2: [kr21], Ea2: [Ea21]}, IVP)

cA1_f.set("cA1_f")
cB1_f.set("cB1_f")
cC1_f.set("cC1_f")
DAG.output([cA1_f, cB1_f, cC1_f])

print(cA1_f, " = ", cA1_f.str())

# %%
# reactor #2 initial concentrations
cC2_0 = DAG.add_var("cC2_0")

# reactor #2 parameters
tau2 = DAG.add_var("tau2")
T2   = DAG.add_var("T2")

# reactor #2 constants
kr12 = DAG.add_var("kr12")
Ea12 = DAG.add_var("Ea12")
kr22 = DAG.add_var("kr22")
Ea22 = DAG.add_var("Ea22")

[cA2_f, cB2_f, cC2_f] = OpIVP({tau: [tau2], T: [T2], cA0: [cA1_f], cB0: [cB1_f], cC0: [cC2_0]},
                              {kr1: [kr12], Ea1: [Ea12], kr2: [kr22], Ea2: [Ea22]}, IVP)

cA2_f.set("cA2_f")
cB2_f.set("cB2_f")
cC2_f.set("cC2_f")
DAG.output([cA2_f, cB2_f, cC2_f])

# %% [markdown]
# Lastly, we define the two process outputs of interest:

# %%
xC1_f = cC1_f / (cA1_f + cB1_f + cC1_f)
xC1_f.set("xC1_f")

print(xC1_f, " = ", xC1_f.str())

# %%
xB2_f = cB2_f / (cA2_f + cB2_f + cC2_f)
xB2_f.set("xB2_f")

print(xB2_f, " = ", xB2_f.str())

# %% [markdown]
# This DAG can be evaluated as any other DAG using the `eval` method; the result agrees with the two successive simulations above:

# %%
inputs = [tau1, T1, tau2, T2, kr11, Ea11, kr21, Ea21, cA1_0, cB1_0, cC1_0, kr12, Ea12, kr22, Ea22, cC2_0]
nominal = [400, 800, 400, 800, 6.66e-3, 2.52e3, 1.03, 5e3, 2, 0, 0, 6.66e-3, 2.52e3, 1.03, 5e3, 0]
[xC1_f_val, xB2_f_val] = DAG.eval([xC1_f, xB2_f], inputs, nominal)
print("xC1_f =", xC1_f_val)
print("xB2_f =", xB2_f_val)
print("from the simulations above: xC1_f =", traj1[-1, 3] / traj1[-1, 1:].sum(), " xB2_f =", traj2[-1, 2] / traj2[-1, 1:].sum())

# %% [markdown]
# Finally, we conduct multiple DAG evaluations over the operation domain $(\tau_j,T_j)\in [250,800] \times [250, 1000]$:

# %%
from scipy.stats import qmc

sampler = qmc.Sobol(d=4, scramble=False, seed=1)
samInput = sampler.random_base2(m=12)  # 2^12 ~ 4,096 evaluations
samInput[:, 0] = 250 + samInput[:, 0] * (800 - 250)
samInput[:, 1] = 250 + samInput[:, 1] * (1000 - 250)
samInput[:, 2] = 250 + samInput[:, 2] * (800 - 250)
samInput[:, 3] = 250 + samInput[:, 3] * (1000 - 250)

print(samInput)

# %%
DAG.options.MAXTHREAD = 0  # use all available threads
samOutput = DAG.veval([xC1_f, xB2_f], [tau1, T1, tau2, T2], samInput.tolist(),
                      [kr11, Ea11, kr21, Ea21, cA1_0, cB1_0, cC1_0, kr12, Ea12, kr22, Ea22, cC2_0],
                      [6.66e-3, 2.52e3, 1.03, 5e3, 2, 0, 0, 6.66e-3, 2.52e3, 1.03, 5e3, 0],
                      walltime=True)

print(np.array(samOutput))

# %%
samOutput = np.array(samOutput)

fig = plt.figure(figsize=(14, 9))
ax = plt.axes(projection='3d')
ax.set_xlabel(r'$\tau_1$')
ax.set_ylabel(r'$T_1$')
ax.set_zlabel(r'$x_{C,1}(\tau_1)$')

ax.scatter(samInput[:, 0], samInput[:, 1], samOutput[:, 0], s=1)
ax.view_init(elev=25, azim=200)

plt.show()

# %% [markdown]
# We can also evaluate the derivatives of the two outputs $x_{C,1}, x_{B,2}$ with respect to any of the inputs $\tau_1,T_1,\tau_2,T_2$. The derivatives of the `FFODESLV` operations are computed by forward or adjoint sensitivity analysis (`FFODESLV.options.GRADIENT`):

# %%
[Row, Col, Grad] = DAG.fdiff([xC1_f, xB2_f], [tau1, T1, tau2, T2])
print("Row:", Row)
print("Col:", Col)

Grad_val = DAG.eval(Grad, inputs, nominal)
print("Grad =", Grad_val)
