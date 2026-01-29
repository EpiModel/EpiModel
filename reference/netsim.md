# Stochastic Network Models

Simulates stochastic network epidemic models for infectious disease.

## Usage

``` r
netsim(x, param, init, control)
```

## Arguments

- x:

  If `control$start == 1`, either a fitted network model object of class
  `netest` or a list of such objects. If `control$start > 1`, an object
  of class `netsim`. When multiple networks are used, the node sets
  (including network size and nodal attributes) are assumed to be the
  same for all networks.

- param:

  Model parameters, as an object of class `param.net`.

- init:

  Initial conditions, as an object of class `init.net`.

- control:

  Control settings, as an object of class `control.net`.

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

## References

Jenness SM, Goodreau SM and Morris M. EpiModel: An R Package for
Mathematical Modeling of Infectious Disease over Networks. Journal of
Statistical Software. 2018; 84(8): 1-47.

## See also

Extract the model results with
[`as.data.frame.netsim()`](http://epimodel.github.io/EpiModel/reference/as.data.frame.icm.md).
Summarize the time-specific model results with
[`summary.netsim()`](http://epimodel.github.io/EpiModel/reference/summary.netsim.md).
Plot the model results with
[`plot.netsim()`](http://epimodel.github.io/EpiModel/reference/plot.netsim.md).

## Examples

``` r
if (FALSE) { # \dontrun{
## Example 1: SI Model without Network Feedback
# Network model estimation
nw <- network_initialize(n = 100)
nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

# Epidemic model
param <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.15)
init <- init.net(i.num = 10, i.num.g2 = 10)
control <- control.net(type = "SI", nsteps = 100, nsims = 5, verbose.int = 0)
mod1 <- netsim(est1, param, init, control)

# Print, plot, and summarize the results
mod1
plot(mod1)
summary(mod1, at = 50)

## Example 2: SIR Model with Network Feedback
# Recalculate dissolution coefficient with departure rate
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20,
                               d.rate = 0.0021)

# Reestimate the model with new coefficient
est2 <- netest(nw, formation, target.stats, coef.diss)

# Reset parameters to include demographic rates
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
} # }
```
