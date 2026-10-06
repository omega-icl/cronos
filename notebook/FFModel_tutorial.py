# %% [markdown]
# # FFModel: analysing a model before solving it
#
# `FFModel` holds the declaration of a model -- its domains, states, inputs, equations and outputs -- and, through
# `setup()`, works out what the model *is* before any solver touches it.  Setup produces a **working model**, in
# which
#
# * derivatives of order higher than one are **order-reduced** to first order, through auxiliary states and `LINK`
#   rows;
# * high-index differential-algebraic equations are **index-reduced** to index one, by differentiating the
#   constraints that need it;
#
# and it **analyses** the result:
#
# * the **classification** of each block of equations -- algebraic, ordinary differential, differential-algebraic,
#   parabolic, elliptic, hyperbolic -- from its principal symbol;
# * the **well-posedness** checks, which flag a model that cannot be solved as declared;
# * the **initial data** a high-index model admits: which initial values may be chosen, and which are fixed by
#   hidden constraints;
# * the **degrees of freedom**: whether the collocated rows balance the unknowns.
#
# The solvers -- `ODESLV` for ODE/DAE models and `OCFESLV` for PDE/DAE models -- are built on `FFModel`, so the same
# declaration and the same analyses come with them.  This tutorial uses `FFModel` on its own, on eight models of
# increasing complexity, and reads each analysis from the model's **report**.

# %%
from pymcpp import FFGraph, FFPartial, FFEval, FFIntegral
from cronos import FFModel, FFDom

OpP, OpE, OpI = FFPartial(), FFEval(), FFIntegral()

Role = FFModel.EqnRole
def role(r): return FFModel.EqnOptions(r)

# %% [markdown]
# ## 1. A system of algebraic equations
#
# The intersection of the unit circle with the line $y = x$: two states without domains, two equations.
#
# A model without domains is **lumped**; its equations are collocated once.

# %%
G = FFGraph()
x, y = G.add_var("x"), G.add_var("y")

M = FFModel(G)
M.add_state(x, [], ref=0.5)
M.add_state(y, [], ref=0.5)
M.add_equation(x * x + y * y - 1.)
M.add_equation(x - y)

print("setup:", M.setup(), "|", M.setup_status)
M.report()

# %% [markdown]
# The report lists the working model -- here identical to the declaration -- followed by the analyses.  The block
# is classified `ALGEBRAIC_LUMPED`, and two rows for two unknowns balance.
#
# Adding a third equation, $x + y = 1.4$, over-specifies the system.  `setup()` still succeeds -- an over-specified
# system can be solved in the least-squares sense -- but both the well-posedness check and the degrees of freedom
# say so:

# %%
M.add_equation(x + y - 1.4)

M.setup()
M.report()

# %% [markdown]
# ## 2. A system of ODEs
#
# An isothermal continuous stirred-tank reactor (CSTR) with the reaction $A \to B$, residence time $\tau = 2$ and
# rate constant $k = 0.5$:
#
# $$\begin{aligned}
# \frac{dC_A}{dt} &= \frac{C_{A,\rm in} - C_A}{\tau} - k\,C_A\\ \frac{dC_B}{dt} &= -\frac{C_B}{\tau} + k\,C_A
# \end{aligned}$$
#
# The evolution domain $t \in [0, 5]$ is discretised into 4 finite elements with 4 Radau (LGR) nodes each.  The
# differential equations hold everywhere but at $t = 0$ (the mask `FFDom.ALL - FFDom.LB`), where the initial conditions hold
# instead (the mask `FFDom.LB`).

# %%
G = FFGraph()
t, CA, CB = G.add_var("t"), G.add_var("CA"), G.add_var("CB")

M = FFModel(G)
M.add_domain(t, FFDom(0., 5., 4, FFDom.LGR, 4))
M.set_evolution_domain(t)
M.add_state(CA, [t], ref=1.)
M.add_state(CB, [t], ref=0.)
M.add_equation(OpP(CA, t) - (1. - CA) / 2. + 0.5 * CA, [t], [FFDom.ALL - FFDom.LB], role(Role.INTERIOR))
M.add_equation(OpP(CB, t) + CB / 2. - 0.5 * CA,        [t], [FFDom.ALL - FFDom.LB], role(Role.INTERIOR))
M.add_equation(CA - 1.,                                [t], [FFDom.LB],             role(Role.INITIAL))
M.add_equation(CB,                                     [t], [FFDom.LB],             role(Role.INITIAL))
#M.add_output(CB, [t], point=[5.])
M.add_output(OpE(CB, t, 5.))

M.setup()
M.report()

# %% [markdown]
# `DIFFERENTIAL_ORDINARY`: every state is differentiated, and the coefficient of the derivatives (here the
# identity) is regular -- index 0.  The output $C_B(5)$ is listed as declared.

# %% [markdown]
# ## 3. An ODE with a state transition
#
# A drug dose: $x' = -k\,x$ with $x(0) = 1$, and a dose $d$ given at $t = 0.4$, so that
# $x(0.4^+) = x(0.4^-) + d$.  `add_transition(left, right, t, tau)` declares that `left` evaluated just before
# $\tau$ equals `right` just after it.

# %%
G = FFGraph()
t, x, d = G.add_var("t"), G.add_var("x"), G.add_var("d")

M = FFModel(G)
M.add_domain(t, FFDom(0., 1., 5, FFDom.LGR, 4))
M.set_evolution_domain(t)
M.add_state(x, [t], ref=1.)
M.add_input(d, ref=0.7)
M.add_equation(OpP(x, t) + 0.8 * x, [t], [FFDom.ALL - FFDom.LB], role(Role.INTERIOR))
M.add_equation(x - 1.,              [t], [FFDom.LB],             role(Role.INITIAL))
M.add_transition(x + d, x, t, 0.4)

M.setup()
M.report()

# %% [markdown]
# The report lists the transition, which setup has validated, and the classification is unchanged: a transition
# restarts the ODE at $\tau$, it does not change its character.

# %% [markdown]
# ## 4. A DAE of index 1
#
# The CSTR again, with a second-order reaction whose rate $r$ is an algebraic state:
#
# $$\begin{aligned}
# \frac{dC_A}{dt} &= -r\\ 0 &= r - k\,C_A^2
# \end{aligned}$$
#
# The algebraic equation holds at every node, $t = 0$ included (the mask `FFDom.ALL`): it determines $r$ from $C_A$
# everywhere, so only $C_A$ needs an initial condition.

# %%
G = FFGraph()
t, CA, r = G.add_var("t"), G.add_var("CA"), G.add_var("r")

M = FFModel(G)
M.add_domain(t, FFDom(0., 5., 4, FFDom.LGR, 4))
M.set_evolution_domain(t)
M.add_state(CA, [t], ref=1.)
M.add_state(r, [t], ref=0.5)
M.add_equation(OpP(CA, t) + r,     [t], [FFDom.ALL - FFDom.LB], role(Role.INTERIOR))
M.add_equation(r - 0.5 * CA * CA,  [t], [FFDom.ALL],            role(Role.INTERIOR))
M.add_equation(CA - 1.,            [t], [FFDom.LB],             role(Role.INITIAL))

M.setup()
M.report()

# %% [markdown]
# `DIFFERENTIAL_ALGEBRAIC` of index 1: the algebraic equation can be solved for $r$ directly, so no reduction is
# needed and the model balances with a single initial condition.

# %% [markdown]
# ## 5. High-index DAEs
#
# ### The pendulum (index 3)
#
# A pendulum of unit length in Cartesian coordinates: positions $(x, y)$, velocities $(u, v)$, and the tension
# $\lambda$ as the Lagrange multiplier of the length constraint:
#
# $$\begin{aligned}
# x' &= u\\
# y' &= v\\
# u' &= -\lambda x\\
# v' &= -\lambda y - g\\
# 0  &= x^2 + y^2 - 1
# \end{aligned}$$
#
# The constraint does not involve $\lambda$: it must be differentiated **twice** before it determines $\lambda$, so
# the system has differential index 3.  We release the pendulum from $x = 1$ with no vertical velocity.

# %%
G = FFGraph()
t = G.add_var("t")
x, y, u, v, lam = [G.add_var(n) for n in ("x", "y", "u", "v", "lam")]

M = FFModel(G)
M.add_domain(t, FFDom(0., 1., 4, FFDom.LGR, 4))
M.set_evolution_domain(t)
for s, s0 in ((x, 1.), (y, 0.), (u, 0.), (v, 0.), (lam, 0.)):
    M.add_state(s, [t], ref=s0)
M.add_equation(OpP(x, t) - u,              [t], [FFDom.ALL - FFDom.LB], role(Role.INTERIOR))
M.add_equation(OpP(y, t) - v,              [t], [FFDom.ALL - FFDom.LB], role(Role.INTERIOR))
M.add_equation(OpP(u, t) + lam * x,        [t], [FFDom.ALL - FFDom.LB], role(Role.INTERIOR))
M.add_equation(OpP(v, t) + lam * y + 9.81, [t], [FFDom.ALL - FFDom.LB], role(Role.INTERIOR))
M.add_equation(x * x + y * y - 1.,         [t], [FFDom.ALL],            role(Role.INTERIOR))
M.add_equation(x - 1.,                     [t], [FFDom.LB],             role(Role.INITIAL))
M.add_equation(v,                          [t], [FFDom.LB],             role(Role.INITIAL))

M.setup()
M.report()

# %% [markdown]
# The working model has replaced the length constraint by its **second derivative**, which involves $\lambda$
# (the method of dummy derivatives): the `INDEX REDUCTION` section says which constraint was differentiated, how
# many times, and which variable it now determines.  The reduced system is of index 1.
#
# The two derivatives of the constraint are not free, though: the constraint itself, $x^2 + y^2 = 1$, and its
# first derivative, $x u + y v = 0$, must still hold at $t = 0$.  These are **hidden constraints**.  The
# `INITIAL DATA` section counts them: of the 4 initial values of $x, y, u, v$, only 2 may be chosen -- our
# $x(0) = 1$ and $v(0) = 0$ -- and the other 2 follow from the hidden constraints.  
#
# The default option `REDUCE.HIDDEN_IC=True` makes setup add the hidden constraints as `INITIAL` rows, and the model balances with
# two initial conditions.

# %%
M.options.REDUCE.HIDDEN_IC = False

M.setup()
M.report()

# %% [markdown]
# The opposite mistake -- prescribing **all four** initial values, as one would for an ODE -- is caught too: with
# the hidden constraints among the rows, two of the four are redundant with them, and the initial point is
# over-determined -- a surplus of 2 rows, which a solver would satisfy only in the least-squares sense.

# %%
M.add_equation(y, [t], [FFDom.LB], role(Role.INITIAL))
M.add_equation(u, [t], [FFDom.LB], role(Role.INITIAL))
M.options.REDUCE.HIDDEN_IC = True

M.setup()
M.report()

# %% [markdown]
# ### A CSTR with a fast equilibrium (index 2)
#
# The reaction $A \rightleftharpoons B$ is fast enough to be at equilibrium, $C_B = K\,C_A$: its rate $r$ is no
# longer given by a kinetic law, but becomes whatever keeps the equilibrium:
#
# $$\begin{aligned}
# \frac{dC_A}{dt} &= \frac{C_{A,\rm in} - C_A}{\tau} - r\\
# \frac{dC_B}{dt} &= -\frac{C_B}{\tau} + r\\
# 0 &= C_B - K\,C_A
# \end{aligned}$$
#
# The equilibrium does not involve $r$; differentiated once, it does -- an index-2 system, common in chemical
# engineering whenever an equilibrium is imposed rather than modelled.

# %%
G = FFGraph()
t, CA, CB, r = G.add_var("t"), G.add_var("CA"), G.add_var("CB"), G.add_var("r")

M = FFModel(G)
M.add_domain(t, FFDom(0., 5., 4, FFDom.LGR, 4))
M.set_evolution_domain(t)
M.add_state(CA, [t], ref=0.5)
M.add_state(CB, [t], ref=1.)
M.add_state(r, [t], ref=0.)
M.add_equation(OpP(CA, t) - (1. - CA) / 2. + r, [t], [FFDom.ALL - FFDom.LB], role(Role.INTERIOR))
M.add_equation(OpP(CB, t) + CB / 2. - r,        [t], [FFDom.ALL - FFDom.LB], role(Role.INTERIOR))
M.add_equation(CB - 2. * CA,                    [t], [FFDom.ALL],            role(Role.INTERIOR))
M.add_equation(CA + CB - 1.5,                   [t], [FFDom.LB],             role(Role.INITIAL))

M.setup()
M.report()

# %% [markdown]
# One differentiation, one hidden constraint -- the equilibrium itself at $t = 0$ -- and so a single initial
# value to choose: here the total concentration $C_A + C_B = 1.5$; the equilibrium then splits it.

# %% [markdown]
# ## 6. A parabolic PDE
#
# Heat conduction in a rod, $u_t = \alpha\,u_{zz}$ on $z \in (0, 1)$, with the temperature fixed at $z = 0$, a heat
# flux $q$ entering at $z = 1$, and a uniform initial temperature.  The outputs are the temperature at the heated
# end at the final time, $u(1,1)$ and the mean temperature, $\iint_{[0,1]^2} u(t,z)\, dt\, dz$.

# %%
G = FFGraph()
t, z, u = G.add_var("t"), G.add_var("z"), G.add_var("u")

M = FFModel(G)
M.add_domain(t, FFDom(0., 1., 4, FFDom.LGR, 3))
M.add_domain(z, FFDom(0., 1., 4, FFDom.LGL, 5))
M.set_evolution_domain(t)
M.add_state(u, [t, z], ref=0.)
M.add_equation(OpP(u, t) - 0.1 * OpP(u, {z: 2}), [t, z], [FFDom.ALL - FFDom.LB, FFDom.ALL - FFDom.LB - FFDom.UB], role(Role.INTERIOR))
M.add_equation(u,                                [t, z], [FFDom.ALL - FFDom.LB, FFDom.LB],                        role(Role.BOUNDARY))
M.add_equation(0.1 * OpP(u, z) - 1.,             [t, z], [FFDom.ALL - FFDom.LB, FFDom.UB],                        role(Role.BOUNDARY))
M.add_equation(u,                                [t, z], [FFDom.LB,             FFDom.ALL - FFDom.LB - FFDom.UB], role(Role.INITIAL))
M.add_output(OpE(u, {t: 1, z: 1}, {t: 1., z: 1.}))    # u(1, 1): the directions, then the point
M.add_output(OpI(u, {t: 1, z: 1}))                    # the space-time mean of u

M.setup()
M.report()

# %% [markdown]
# **Order reduction**: the second derivative $u_{zz}$ became the first derivative of an auxiliary state
# $D_z u$, defined by the `LINK` row $\partial_z u - D_z u = 0$.  The flux condition at $z = 1$ is expressed with it
# too.  **Classification**: `PARABOLIC` -- the evolution coefficient is singular (the `LINK` row has no time
# derivative), which is the descriptor form a parabolic equation takes once order-reduced.  The `STRUCTURAL INDEX`
# section gives the index per direction.  The two outputs are listed as declared; how a solver evaluates them is its
# own business.

# %% [markdown]
# ## 7. An elliptic PDE (no evolution domain)
#
# The Poisson equation $u_{xx} + u_{yy} + 1 = 0$ on the unit square, with $u = 0$ on all four sides.  There is no
# evolution domain: the problem is a boundary-value problem in both directions.

# %%
G = FFGraph()
x, y, u = G.add_var("x"), G.add_var("y"), G.add_var("u")

M = FFModel(G)
M.add_domain(x, FFDom(0., 1., 2, FFDom.LGL, 5))
M.add_domain(y, FFDom(0., 1., 2, FFDom.LGL, 5))
M.add_state(u, [x, y], ref=0.)
M.add_equation(OpP(u, {x: 2}) + OpP(u, {y: 2}) + 1., [x, y], [FFDom.ALL - FFDom.LB - FFDom.UB, FFDom.ALL - FFDom.LB - FFDom.UB], role(Role.INTERIOR))
M.add_equation(u,                                    [x, y], [FFDom.ALL - FFDom.LB - FFDom.UB, FFDom.LB],                        role(Role.BOUNDARY))
M.add_equation(u,                                    [x, y], [FFDom.ALL - FFDom.LB - FFDom.UB, FFDom.UB],                        role(Role.BOUNDARY))
M.add_equation(u,                                    [x, y], [FFDom.LB,                        FFDom.ALL],                       role(Role.BOUNDARY))
M.add_equation(u,                                    [x, y], [FFDom.UB,                        FFDom.ALL],                       role(Role.BOUNDARY))

M.setup()
M.report()

# %% [markdown]
# Both second derivatives are order-reduced, with one auxiliary and one `LINK` row each, and the block is
# classified `ELLIPTIC`.  For an order-reduced system this takes some care: the symbol of the first derivatives
# alone is singular in every direction (the two `LINK` rows are parallel), and it is the *weighted* symbol, in which
# the `LINK` rows' own auxiliaries count, that is nonsingular -- its determinant is $\xi_x^2 + \xi_y^2$.

# %% [markdown]
# ## 8. A hyperbolic PDE
#
# Linear advection $u_t + a\,u_z = 0$ with speed $a = 1.5 > 0$: information travels from $z = 0$ towards $z = 1$, so
# the boundary condition belongs at the **inflow** boundary $z = 0$, and the equation itself holds up to and
# including the outflow boundary $z = 1$.

# %%
G = FFGraph()
t, z, u = G.add_var("t"), G.add_var("z"), G.add_var("u")

M = FFModel(G)
M.add_domain(t, FFDom(0., 1., 4, FFDom.LGR, 3))
M.add_domain(z, FFDom(0., 1., 4, FFDom.LGL, 4))
M.set_evolution_domain(t)
M.add_state(u, [t, z], ref=0.)
M.add_equation(OpP(u, t) + 1.5 * OpP(u, z), [t, z], [FFDom.ALL - FFDom.LB, FFDom.ALL - FFDom.LB], role(Role.INTERIOR))
M.add_equation(u - 1.,                      [t, z], [FFDom.ALL - FFDom.LB, FFDom.LB],             role(Role.BOUNDARY))
M.add_equation(u,                           [t, z], [FFDom.LB,             FFDom.ALL],            role(Role.INITIAL))

M.setup()
M.report()

# %% [markdown]
# `EVOL_HYPERBOLIC`: the problem evolves in $t$ along real characteristics in $z$, and the condition at the inflow
# boundary is where the characteristics require it.  Moved to the **outflow** side, the condition no longer makes
# sense: information leaves the domain there, and the inflow is left unspecified.  The well-posedness analysis
# checks the boundary conditions against the characteristics, face by face:

# %%
G = FFGraph()
t, z, u = G.add_var("t"), G.add_var("z"), G.add_var("u")

M = FFModel(G)
M.add_domain(t, FFDom(0., 1., 4, FFDom.LGR, 3))
M.add_domain(z, FFDom(0., 1., 4, FFDom.LGL, 4))
M.set_evolution_domain(t)
M.add_state(u, [t, z], ref=0.)
M.add_equation(OpP(u, t) + 1.5 * OpP(u, z), [t, z], [FFDom.ALL - FFDom.LB, FFDom.ALL - FFDom.LB], role(Role.INTERIOR))
M.add_equation(u - 1.,                      [t, z], [FFDom.ALL - FFDom.LB, FFDom.UB],             role(Role.BOUNDARY))
M.add_equation(u,                           [t, z], [FFDom.LB,             FFDom.ALL],            role(Role.INITIAL))

M.setup()
M.report()

# %% [markdown]
# Two findings: the condition sits at a face where **no characteristic enters** -- data at the wrong end -- and
# the rows no longer balance the unknowns.

# %% [markdown]
# ## Summary
#
# | model | classification | what the analyses showed |
# |---|---|---|
# | circle ∩ line | `ALGEBRAIC_LUMPED` | an extra equation: over-specified, flagged |
# | CSTR | `DIFFERENTIAL_ORDINARY` | regular, balanced |
# | dose | `DIFFERENTIAL_ORDINARY` | the transition, validated at setup |
# | CSTR with algebraic rate | `DIFFERENTIAL_ALGEBRAIC` | index 1: no reduction, one initial condition |
# | pendulum | `DIFFERENTIAL_ALGEBRAIC` | index 3 reduced; 2 hidden constraints; too few or too many initial values flagged |
# | CSTR at equilibrium | `DIFFERENTIAL_ALGEBRAIC` | index 2 reduced; 1 hidden constraint |
# | heated rod | `PARABOLIC` | order reduction of $u_{zz}$ |
# | Poisson | `ELLIPTIC` | order reduction in both directions; a missing side flagged |
# | advection | `EVOL_HYPERBOLIC` | a boundary condition on the outflow side flagged |
#
# A model declared this way is passed to a solver as it is: `ODESLV` for the lumped and ODE/DAE models, `OCFESLV`
# for all of them, including the PDEs.  See the `ODESLV` and `OCFESLV` tutorials.
