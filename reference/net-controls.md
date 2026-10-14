# Control Settings for Network Model Estimation and Simulation

[`netest()`](https://epimodel.github.io/EpiModel/reference/netest.md),
[`netdx()`](https://epimodel.github.io/EpiModel/reference/netdx.md), and
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md)
pass control settings to the `ergm` and `tergm` packages, which estimate
and simulate the network model, through their `set.control.ergm`,
`set.control.tergm`, and `set.control.ergm.ego` arguments. This page
lists the control function that each argument takes, the defaults that
EpiModel sets, how to use
[`netdx()`](https://epimodel.github.io/EpiModel/reference/netdx.md) to
find the stage that needs adjusting, and the settings most often
changed.

The worked examples are in the [Improving Model Estimation and
Simulation](https://epimodel.github.io/sismid/11_advanced/mod11-NetFitting.html)
appendix of the Network Modeling for Epidemics course materials, which
covers full STERGM estimation for short partnership durations, the
estimation settings for large and slow ERGM fits, and the simulation
settings for diagnostics and epidemic models.

## Where Each Setting Applies

Each argument takes an object built by a specific control function.
Passing an object from another one stops with an error that names the
right function.

- **Estimation in
  [`netest()`](https://epimodel.github.io/EpiModel/reference/netest.md):**

  - `set.control.ergm` takes
    [`ergm::control.ergm()`](https://rdrr.io/pkg/ergm/man/control.ergm.html),
    for the ERGM fit with the edges dissolution approximation
    (`edapprox = TRUE`).

  - `set.control.ergm.ego` takes
    [`ergm.ego::control.ergm.ego()`](https://rdrr.io/pkg/ergm.ego/man/control.ergm.ego.html),
    for the same fit when `nw` is an `egor` object.

  - `set.control.tergm` takes
    [`tergm::control.tergm()`](https://rdrr.io/pkg/tergm/man/control.tergm.html),
    for the full STERGM fit (`edapprox = FALSE`).

- **Simulation in
  [`netdx()`](https://epimodel.github.io/EpiModel/reference/netdx.md)
  and
  [`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md):**

  - `set.control.ergm` takes
    [`ergm::control.simulate.formula()`](https://rdrr.io/pkg/ergm/man/control.simulate.ergm.html),
    for the cross-sectional network that the dynamic simulation starts
    from, and for the static diagnostics of `netdx(dynamic = FALSE)`. It
    is used only for fits with `edapprox = TRUE`; a full STERGM fit
    starts from the network saved by
    [`netest()`](https://epimodel.github.io/EpiModel/reference/netest.md).

  - `set.control.tergm` takes
    [`tergm::control.simulate.formula.tergm()`](https://rdrr.io/pkg/tergm/man/control.simulate.tergm.html),
    for each time step of the dynamic network simulation, both in the
    diagnostics and in the epidemic model.

In a model with several network layers,
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md)
takes one control object per layer through
[`multilayer()`](https://epimodel.github.io/EpiModel/reference/multilayer.md).

## Defaults Set by EpiModel

Three defaults differ from those of the control functions themselves:

- `set.control.tergm` in
  [`netest()`](https://epimodel.github.io/EpiModel/reference/netest.md)
  is `control.tergm(MCMC.maxchanges = .Machine$integer.max)`.

- `set.control.tergm` in
  [`netdx()`](https://epimodel.github.io/EpiModel/reference/netdx.md)
  and
  [`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md)
  is
  `control.simulate.formula.tergm(MCMC.maxchanges = .Machine$integer.max)`.
  This raises the number of network changes that `tergm` can record in
  one simulation above its default limit of 1e6, which a long
  [`netdx()`](https://epimodel.github.io/EpiModel/reference/netdx.md)
  run can exceed.

- `set.control.ergm` in
  [`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md)
  is `control.simulate.formula(MCMC.burnin = 2e5)`, a longer burn-in for
  the starting network than the `ergm` default.

A control object passed to one of these arguments replaces the default
as a whole, not only the settings it names, so restate these settings
alongside any others (see the examples).

## Which Stage to Adjust

An error or a stalled fit in
[`netest()`](https://epimodel.github.io/EpiModel/reference/netest.md)
calls for the estimation settings below. When the fit succeeds but the
dynamic diagnostics of
[`netdx()`](https://epimodel.github.io/EpiModel/reference/netdx.md) miss
the formation targets, the cause can lie in the simulation, in the edges
dissolution approximation, or in the cross-sectional fit:

1.  **The per-step burn-in of the simulation.** By default, `tergm` sets
    the number of proposals in each time step adaptively. With short
    partnership durations it can end a time step before the network has
    replaced the partnerships that dissolved, so statistics such as the
    edge count settle below their targets. This applies to fits with and
    without the approximation. Setting `MCMC.burnin.min` and
    `MCMC.burnin.max` in `set.control.tergm` to the same value (such as
    1e4) fixes the number of proposals per time step and is quick to
    try.

2.  **The edges dissolution approximation.** For fits with
    `edapprox = TRUE`, run static diagnostics with
    `netdx(dynamic = FALSE)`. If they are on target but the dynamic
    diagnostics remain off with a fixed per-step burn-in, the
    approximation is the likely cause. It loses accuracy when many
    partnerships form and dissolve at each time step, that is, when
    durations are short or density is high. Fit the full STERGM with
    `edapprox = FALSE`, or redefine the time step as a shorter unit so
    that durations span more time steps, rescaling every other
    time-based parameter of the model to match.

3.  **The cross-sectional fit.** If the static diagnostics are also off,
    the ERGM fit has not reached the targets, or the static simulation
    has not moved far enough from its starting network. On a large
    network, raise `MCMC.burnin` and `MCMC.interval` in the
    `set.control.ergm` of
    [`netdx()`](https://epimodel.github.io/EpiModel/reference/netdx.md)
    before changing the estimation settings of
    [`netest()`](https://epimodel.github.io/EpiModel/reference/netest.md).

Read the `Pct Diff` column of
[`print.netdx()`](https://epimodel.github.io/EpiModel/reference/print.netdx.md)
alongside the `Z Score`. The standard error behind the Z score shrinks
with the effective sample size, so many static draws can make a
difference of a fraction of a percent statistically significant, while
the autocorrelation between time steps keeps the effective sample size
of a dynamic simulation, and so its Z scores, smaller.

The simulation settings that bring the
[`netdx()`](https://epimodel.github.io/EpiModel/reference/netdx.md)
diagnostics on target belong in
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md)
as well: the epidemic simulation steps the network forward with the same
`tergm` simulation, so without them it runs on the network that the
diagnostics flagged.

## Settings Most Often Changed

The control functions have many arguments, of which a few account for
most adjustments in practice. Their help pages document the rest.

For the ERGM fit, in
[`ergm::control.ergm()`](https://rdrr.io/pkg/ergm/man/control.ergm.html):

- `MCMLE.maxit`: the maximum number of MCMLE iterations. Raise it for a
  fit that is still moving toward the targets when it stops.

- `MCMC.samplesize`, `MCMC.interval`, and `MCMC.burnin`: the size,
  spacing, and burn-in of the MCMC sample drawn at each iteration.
  Larger networks generally need a longer interval and burn-in. By
  default `ergm` extends the sample until it reaches the effective size
  set by `MCMLE.effectiveSize`; setting `MCMLE.effectiveSize = NULL`
  holds it at `MCMC.samplesize` draws, which bounds the time per
  iteration.

- `SAN.nsteps.times` and `SAN.maxit`: simulated annealing (SAN) builds
  the starting network that matches the target statistics.
  `SAN.nsteps.times` scales the total number of SAN proposals, and
  `SAN.maxit` divides them into that many runs. Raise `SAN.nsteps.times`
  when SAN does not reach the targets.

- `MCMC.prop`: hints for the MCMC proposal (see
  [ergm::ergmHint](https://rdrr.io/pkg/ergm/man/ergmHint.html)). When
  most ties fall within a few cells of a mixing matrix, a stratified
  proposal such as `~strat(attr = ~age.grp) + sparse` stops the sampler
  from spending its proposals on cells that hold few or no ties. Listing
  any hint replaces the default ones, so list `sparse` as well.

- `MCMC.maxedges`: the maximum number of edges in a sampled network,
  unlimited by default. A finite value makes a fit that wanders into
  dense, degenerate networks stop with an error rather than run for
  hours.

- `MCMLE.steplength` with `MCMLE.steplength.margin = NULL`: a fixed step
  length below 1 (such as 0.5) keeps the first iterations of a hard fit
  from overshooting. With the default margin, `MCMLE.steplength` is only
  the maximum step length.

- `parallel`: the number of cores across which the MCMC sample is drawn.

For the full STERGM fit, in
[`tergm::control.tergm()`](https://rdrr.io/pkg/tergm/man/control.tergm.html):

- `EGMME.MCMC.burnin.min` and `EGMME.MCMC.burnin.max`: the number of
  proposals per time step during estimation, set adaptively by default.
  Setting both to the same value fixes it, which can stabilize a fit
  that otherwise fails. In the course examples, 1e4 works for a network
  of 1,000 nodes, and larger networks may need a value on the order of
  the network size. Fix the per-step burn-in of the simulation to the
  same value.

- `SAN.nsteps.times`, `SAN.maxit`, and `MCMC.prop` act as in
  [`ergm::control.ergm()`](https://rdrr.io/pkg/ergm/man/control.ergm.html).

For the dynamic simulation, in
[`tergm::control.simulate.formula.tergm()`](https://rdrr.io/pkg/tergm/man/control.simulate.tergm.html):

- `MCMC.burnin.min` and `MCMC.burnin.max`: the number of proposals per
  time step. Setting both to the same value fixes it, as described
  above. `tergm` has no `MCMC.interval`; this pair takes its place.

- `MCMC.maxchanges`: restate `.Machine$integer.max` when passing any
  other setting.

- `MCMC.prop`: the default, `~discord + sparse`, proposes a share of its
  toggles (set by `MCMC.prop.args = list(discordance_fraction = )`, 0.5
  by default) on dyads that have already changed within the time step.
  For strongly structured mixing, add a `strat()` hint as in estimation
  and keep `discord` and `sparse` in the formula.

For the cross-sectional simulation, in
[`ergm::control.simulate.formula()`](https://rdrr.io/pkg/ergm/man/control.simulate.ergm.html):

- `MCMC.burnin` and `MCMC.interval`: the burn-in before the starting
  network (or the first static draw) and the spacing between static
  draws.

## References

Krivitsky PN, Hunter DR, Morris M, Klumb C. ergm 4: New Features for
Analyzing Exponential-Family Random Graph Models. Journal of Statistical
Software. 2023; 105(6): 1-44.
[doi:10.18637/jss.v105.i06](https://doi.org/10.18637/jss.v105.i06) .

Krivitsky PN, Handcock MS. A Separable Model for Dynamic Networks.
Journal of the Royal Statistical Society, Series B. 2014; 76(1): 29-46.

## See also

[`netest()`](https://epimodel.github.io/EpiModel/reference/netest.md),
[`netdx()`](https://epimodel.github.io/EpiModel/reference/netdx.md), and
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md),
whose arguments these settings are passed through;
[`print.netdx()`](https://epimodel.github.io/EpiModel/reference/print.netdx.md)
and
[`plot.netdx()`](https://epimodel.github.io/EpiModel/reference/plot.netdx.md)
for reading the diagnostics.

## Examples

``` r
# \donttest{
nw <- network_initialize(n = 100)
formation <- ~edges + concurrent
target.stats <- c(50, 20)
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
est <- netest(nw, formation, target.stats, coef.diss)
#> Starting simulated annealing (SAN)
#> Iteration 1 of at most 4
#> Finished simulated annealing
#> Starting maximum pseudolikelihood estimation (MPLE):
#> Obtaining the responsible dyads.
#> Evaluating the predictor and response matrix.
#> Maximizing the pseudolikelihood.
#> Finished MPLE.
#> Starting Monte Carlo maximum likelihood estimation (MCMLE):
#> Iteration 1 of at most 60:
#> Warning: ‘glpk’ selected as the solver, but package ‘Rglpk’ is not available; falling back to ‘lpSolveAPI’. This should be fine unless the sample size and/or the number of parameters is very big.
#> 1 
#> Optimizing with step length 1.0000.
#> The log-likelihood improved by 1.5495.
#> Estimating equations are not within tolerance region.
#> Iteration 2 of at most 60:
#> 1 
#> Optimizing with step length 1.0000.
#> The log-likelihood improved by 0.0794.
#> Convergence test p-value: 0.0001. 
#> Converged with 99% confidence.
#> Finished MCMLE.
#> This model was fit using MCMC.  To examine model diagnostics and check
#> for degeneracy, use the mcmc.diagnostics() function.

# Fix the number of proposals per time step in the dynamic simulation,
# restating the MCMC.maxchanges default that EpiModel sets
sim.control <- control.simulate.formula.tergm(
  MCMC.burnin.min = 1e4,
  MCMC.burnin.max = 1e4,
  MCMC.maxchanges = .Machine$integer.max
)
dx <- netdx(est, nsims = 2, nsteps = 100, set.control.tergm = sim.control,
            verbose = FALSE)
dx
#> EpiModel Network Diagnostics
#> =======================
#> Diagnostic Method: Dynamic
#> Simulations: 2
#> Time Steps per Sim: 100
#> 
#> Formation Diagnostics
#> ----------------------- 
#>            Target Sim Mean Pct Diff Sim SE Z Score SD(Sim Means) SD(Statistic)
#> edges          50   49.875   -0.250  1.436  -0.087         1.209         5.883
#> concurrent     20   20.885    4.425  1.054   0.840         1.266         4.530
#> 
#> Duration Diagnostics
#> ----------------------- 
#>       Target Sim Mean Pct Diff Sim SE Z Score SD(Sim Means) SD(Statistic)
#> edges     10    9.599   -4.014  0.305  -1.316         0.248         1.152
#> 
#> Dissolution Diagnostics
#> ----------------------- 
#>       Target Sim Mean Pct Diff Sim SE Z Score SD(Sim Means) SD(Statistic)
#> edges    0.1      0.1    0.481  0.003   0.161         0.006         0.042

# Pass the same settings to the epidemic simulation
param <- param.net(inf.prob = 0.3, act.rate = 1)
init <- init.net(i.num = 10)
control <- control.net(type = "SI", nsims = 1, nsteps = 50,
                       set.control.tergm = sim.control, verbose = FALSE)
sim <- netsim(est, param, init, control)
# }

if (FALSE) { # \dontrun{
# Estimation settings for a larger or slower ERGM fit
est <- netest(nw, formation, target.stats, coef.diss,
              set.control.ergm = control.ergm(MCMLE.maxit = 500,
                                              MCMC.samplesize = 1e4,
                                              MCMC.interval = 5e3,
                                              SAN.nsteps.times = 16,
                                              parallel = 4))

# Full STERGM fit with a fixed number of proposals per time step
est <- netest(nw, formation, target.stats, coef.diss, edapprox = FALSE,
              set.control.tergm = control.tergm(
                EGMME.MCMC.burnin.min = 1e4,
                EGMME.MCMC.burnin.max = 1e4,
                MCMC.maxchanges = .Machine$integer.max
              ))
dx <- netdx(est, nsims = 5, nsteps = 500, set.control.tergm = sim.control)
} # }
```
