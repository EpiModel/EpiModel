# Stochastic Network Models

Simulates stochastic network epidemic models for infectious disease.
This is typically the final step in the network modeling pipeline, after
network estimation with
[`netest`](http://epimodel.github.io/EpiModel/reference/netest.md) and
model diagnostics with
[`netdx`](http://epimodel.github.io/EpiModel/reference/netdx.md).

## Usage

``` r
netsim(x, param, init, control)
```

## Arguments

- x:

  If `control$start == 1`, either a fitted network model object of class
  `netest` or a list of such objects. If `control$start > 1`, an object
  of class `netsim`. When multiple networks are used (multi-layer
  models), pass a list of `netest` objects, one per network layer; the
  node sets (including network size and nodal attributes) are assumed to
  be the same for all networks.

- param:

  Model parameters, as an object of class
  [`param.net`](http://epimodel.github.io/EpiModel/reference/param.net.md).
  Includes transmission probability (`inf.prob`), act rate (`act.rate`),
  recovery rate (`rec.rate`), and demographic rates for models with
  vital dynamics. Custom parameters for extended models may also be
  passed through `param.net`.

- init:

  Initial conditions, as an object of class
  [`init.net`](http://epimodel.github.io/EpiModel/reference/init.net.md).
  Specifies the initial number of infected nodes (`i.num`) and, for SIR
  models, recovered nodes (`r.num`). For two-group models, the
  corresponding `.g2` parameters are also required.

- control:

  Control settings, as an object of class
  [`control.net`](http://epimodel.github.io/EpiModel/reference/control.net.md).
  Key settings include `type` (disease model: `"SI"`, `"SIR"`, or
  `"SIS"`), `nsteps` (number of time steps), `nsims` (number of
  simulations), `tergmLite` (lightweight mode for performance), and
  `resimulate.network` (required for models with vital dynamics). For
  extended models, custom module functions are also passed here.

## Value

A list of class `netsim` with the following elements:

- **param:** the epidemic parameters passed into the model through
  `param`, with additional parameters added as necessary.

- **control:** the control settings passed into the model through
  `control`, with additional controls added as necessary.

- **epi:** a list of data frames, one for each epidemiological output
  from the model. Outputs for base models always include the size of
  each compartment, as well as flows in, out of, and between
  compartments.

- **stats:** a list containing two sublists, `nwstats` for any network
  statistics saved in the simulation, and `transmat` for the
  transmission matrix saved in the simulation. See
  [`control.net()`](http://epimodel.github.io/EpiModel/reference/control.net.md)
  for further details.

- **network:** a list of lists of `networkDynamic` or `networkLite`
  objects, with one list of objects for each model simulation.

If `control$raw.output == TRUE`: A list of the raw (pre-processed)
`netsim_dat` objects, for use in simulation continuation.

The `epi` data can be extracted as a data frame with
[`as.data.frame.netsim()`](http://epimodel.github.io/EpiModel/reference/as.data.frame.icm.md),
with options for per-simulation values (`out = "vals"`), means
(`out = "mean"`), standard deviations (`out = "sd"`), or quantiles
(`out = "qnt"`).

## Details

Stochastic network models explicitly represent phenomena within and
across edges (pairs of nodes that remain connected) over time. This
enables edges to have duration, allowing for repeated
transmission-related acts within the same dyad, specification of edge
formation and dissolution rates, control over the temporal sequencing of
multiple edges, and specification of network-level features. A detailed
description of these models, along with examples, is found in the
[Network Modeling for Epidemics](https://epimodel.github.io/sismid/)
course materials.

The `netsim` function performs modeling of both the base model types and
original models. Base model types include one-group and two-group models
with disease types for Susceptible-Infected (SI),
Susceptible-Infected-Recovered (SIR), and
Susceptible-Infected-Susceptible (SIS).

Original models may be parameterized by writing new process modules that
either take the place of existing modules (for example, disease
recovery), or supplement the set of existing processes with a new one
contained in a new module. This functionality is documented in the
[Extending
EpiModel](https://epimodel.github.io/sismid/9_extending/mod9-Intro.html)
section of the [Network Modeling for
Epidemics](https://epimodel.github.io/sismid/) course materials. The
list of modules within `netsim` available for modification is listed in
[`modules.net()`](http://epimodel.github.io/EpiModel/reference/modules.net.md).

## Module Pipeline

At each time step, `netsim` executes a series of modules in sequence.
For base models, the default pipeline is:

1.  **Network resimulation**
    ([`resim_nets`](http://epimodel.github.io/EpiModel/reference/resim_nets.md)):
    updates the network structure (if `resimulate.network = TRUE`).

2.  **Infection**
    ([`infection.net`](http://epimodel.github.io/EpiModel/reference/infection.net.md)):
    simulates disease transmission across discordant edges (where one
    partner is susceptible and the other is infected).

3.  **Recovery**
    ([`recovery.net`](http://epimodel.github.io/EpiModel/reference/recovery.net.md)):
    simulates recovery from infection (SIR and SIS models only).

4.  **Departures**
    ([`departures.net`](http://epimodel.github.io/EpiModel/reference/departures.net.md)):
    simulates node exits from the network (if vital dynamics are
    enabled).

5.  **Arrivals**
    ([`arrivals.net`](http://epimodel.github.io/EpiModel/reference/arrivals.net.md)):
    simulates new node entries into the network (if vital dynamics are
    enabled).

6.  **Prevalence**
    ([`prevalence.net`](http://epimodel.github.io/EpiModel/reference/prevalence.net.md)):
    calculates summary statistics.

See
[`modules.net()`](http://epimodel.github.io/EpiModel/reference/modules.net.md)
for full details on each module.

## Performance and tergmLite

Setting `tergmLite = TRUE` in
[`control.net`](http://epimodel.github.io/EpiModel/reference/control.net.md)
uses a lightweight network representation (`networkLite`) that is
substantially faster than the full `networkDynamic` representation,
often by a factor of 20–50x. This is recommended for large networks or
when running many simulations. The trade-off is that full
`networkDynamic` objects are not available for post-hoc analysis; use
the cumulative edgelist (`cumulative.edgelist = TRUE`) instead to track
partnership histories.

The `resimulate.network` control must be set to `TRUE` when demographic
processes (arrivals and departures) change the network composition over
time. Without it, the network structure evolves independently of the
epidemic and demographic dynamics.

## Multi-Network Models

For models with multiple overlapping network layers (e.g., sexual and
needle-sharing networks), pass a list of
[`netest`](http://epimodel.github.io/EpiModel/reference/netest.md)
objects to the `x` argument, one per network layer. Each layer has its
own formation/dissolution dynamics but shares the same node set. See the
`multilayer` documentation and `test-multinets.R` for examples.

## Restarting and Checkpointing

Simulations can be checkpointed and restarted if interrupted. Set
`.checkpoint.steps` and `.checkpoint.dir` in
[`control.net`](http://epimodel.github.io/EpiModel/reference/control.net.md)
to enable automatic checkpointing. To restart a simulation from a prior
`netsim` output, pass the `netsim` object as `x` and set `control$start`
to one greater than the final time step of the prior simulation. See the
Checkpointing Simulations section of
[`control.net`](http://epimodel.github.io/EpiModel/reference/control.net.md)
for full details.

## References

Jenness SM, Goodreau SM and Morris M. EpiModel: An R Package for
Mathematical Modeling of Infectious Disease over Networks. Journal of
Statistical Software. 2018; 84(8): 1-47.

## See also

Estimate the network model with
[`netest`](http://epimodel.github.io/EpiModel/reference/netest.md) and
diagnose model fit with
[`netdx`](http://epimodel.github.io/EpiModel/reference/netdx.md) before
running simulations. Extract model results with
[`as.data.frame.netsim()`](http://epimodel.github.io/EpiModel/reference/as.data.frame.icm.md).
Summarize the time-specific model results with
[`summary.netsim()`](http://epimodel.github.io/EpiModel/reference/summary.netsim.md).
Plot the model results with
[`plot.netsim()`](http://epimodel.github.io/EpiModel/reference/plot.netsim.md).
Extract the network with
[`get_network`](http://epimodel.github.io/EpiModel/reference/get_network.md),
the transmission matrix with
[`get_transmat`](http://epimodel.github.io/EpiModel/reference/get_transmat.md),
and derive new epi statistics with
[`mutate_epi`](http://epimodel.github.io/EpiModel/reference/mutate_epi.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## Example 1: SI Model without Network Feedback
# Network model estimation
nw <- network_initialize(n = 100)
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

# Epidemic model
param <- param.net(inf.prob = 0.3)
init <- init.net(i.num = 10)
control <- control.net(type = "SI", nsteps = 100, nsims = 5, verbose.int = 0)
mod1 <- netsim(est1, param, init, control)

# Print, plot, and summarize the results
mod1
plot(mod1)
summary(mod1, at = 50)

## Example 2: SIR Model with Network Feedback (Two-Group)
nw <- network_initialize(n = 100)
nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
formation <- ~edges
target.stats <- 50

# Recalculate dissolution coefficient with departure rate
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20,
                               d.rate = 0.0021)

# Reestimate the model with new coefficient
est2 <- netest(nw, formation, target.stats, coef.diss)

# Set parameters to include demographic rates
param <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.15,
                   rec.rate = 0.02, rec.rate.g2 = 0.02,
                   a.rate = 0.002, a.rate.g2 = NA,
                   ds.rate = 0.001, ds.rate.g2 = 0.001,
                   di.rate = 0.001, di.rate.g2 = 0.001,
                   dr.rate = 0.001, dr.rate.g2 = 0.001)
init <- init.net(i.num = 10, i.num.g2 = 10,
                 r.num = 0, r.num.g2 = 0)
control <- control.net(type = "SIR", nsteps = 100, nsims = 5,
                       resimulate.network = TRUE, tergmLite = TRUE)

# Simulate the model with new network fit
mod2 <- netsim(est2, param, init, control)

# Print, plot, and summarize the results
mod2
plot(mod2)
summary(mod2, at = 40)

## Example 3: Post-Simulation Analysis
# Extract epi data as a data frame
as.data.frame(mod1)
as.data.frame(mod1, out = "mean")

# Extract the transmission matrix from simulation 1
get_transmat(mod1, sim = 1)

# Derive new epi statistics and plot them
mod1 <- mutate_epi(mod1, prev = i.num / num)
plot(mod1, y = "prev")
} # }
```
