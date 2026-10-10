# Epidemic Parameters for Stochastic Network Models

Sets the epidemic parameters for stochastic network models simulated
with
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md).

## Usage

``` r
param.net(
  inf.prob,
  inter.eff,
  inter.start,
  act.rate,
  rec.rate,
  a.rate,
  ds.rate,
  di.rate,
  dr.rate,
  inf.prob.g2,
  rec.rate.g2,
  a.rate.g2,
  ds.rate.g2,
  di.rate.g2,
  dr.rate.g2,
  ...
)
```

## Arguments

- inf.prob:

  Probability of infection per transmissible act between a susceptible
  and an infected person. In two-group models, this is the probability
  of infection to the group 1 nodes. This may also vary with the
  infected partner's duration of infection, given as a
  [`by_infection_duration()`](https://epimodel.github.io/EpiModel/reference/by_infection_duration.md)
  object (see Parameters by Duration of Infection below). In multi-layer
  models, this may be a
  [`multilayer()`](https://epimodel.github.io/EpiModel/reference/multilayer.md)
  object with one entry (a probability or a `by_infection_duration`
  object) per network layer (see Multi-Layer Parameters below).

- inter.eff:

  Efficacy of an intervention which affects the per-act probability of
  infection. Efficacy is defined as 1 - the relative hazard of infection
  given exposure to the intervention, compared to no exposure.

- inter.start:

  Time step at which the intervention starts, between 1 and the number
  of time steps specified in the model. This will default to 1 if
  `inter.eff` is defined but this parameter is not.

- act.rate:

  Average number of transmissible acts *per partnership* per unit time
  (see `act.rate` Parameter below). This may also vary with the infected
  partner's duration of infection, given as a
  [`by_infection_duration()`](https://epimodel.github.io/EpiModel/reference/by_infection_duration.md)
  object (see Parameters by Duration of Infection below), or be a
  [`multilayer()`](https://epimodel.github.io/EpiModel/reference/multilayer.md)
  object with one entry per network layer (see Multi-Layer Parameters
  below).

- rec.rate:

  Average rate of recovery with immunity (in `SIR` models) or
  re-susceptibility (in `SIS` models). The recovery rate is the
  reciprocal of the disease duration. For two-group models, this is the
  recovery rate for group 1 persons only. This parameter is only used
  for `SIR` and `SIS` models. This may also vary with the duration of
  infection, given as a
  [`by_infection_duration()`](https://epimodel.github.io/EpiModel/reference/by_infection_duration.md)
  object (see Parameters by Duration of Infection below).

- a.rate:

  Arrival or entry rate. For one-group models, the arrival rate is the
  rate of new arrivals per person per unit time. For two-group models,
  the arrival rate is parameterized as a rate per group 1 person per
  unit time, with the `a.rate.g2` rate set as described below.

- ds.rate:

  Departure or exit rate for susceptible persons. For two-group models,
  it is the rate for group 1 susceptible persons only.

- di.rate:

  Departure or exit rate for infected persons. For two-group models, it
  is the rate for group 1 infected persons only.

- dr.rate:

  Departure or exit rate for recovered persons. For two-group models, it
  is the rate for group 1 recovered persons only. This parameter is only
  used for `SIR` models.

- inf.prob.g2:

  Probability of transmission given a transmissible act between a
  susceptible group 2 person and an infected group 1 person. It is the
  probability of transmission to group 2 members. As for `inf.prob`, it
  may be a
  [`by_infection_duration()`](https://epimodel.github.io/EpiModel/reference/by_infection_duration.md)
  or
  [`multilayer()`](https://epimodel.github.io/EpiModel/reference/multilayer.md)
  object.

- rec.rate.g2:

  Average rate of recovery with immunity (in `SIR` models) or
  re-susceptibility (in `SIS` models) for group 2 persons. This
  parameter is only used for two-group `SIR` and `SIS` models. As for
  `rec.rate`, it may be a
  [`by_infection_duration()`](https://epimodel.github.io/EpiModel/reference/by_infection_duration.md)
  object.

- a.rate.g2:

  Arrival or entry rate for group 2. This may either be specified
  numerically as the rate of new arrivals per group 2 person per unit
  time, or as `NA`, in which case the group 1 rate, `a.rate`, governs
  the group 2 rate. The latter is used when, for example, the first
  group is conceptualized as female, and the female population size
  determines the arrival rate. Such arrivals are evenly allocated
  between the two groups.

- ds.rate.g2:

  Departure or exit rate for group 2 susceptible persons.

- di.rate.g2:

  Departure or exit rate for group 2 infected persons.

- dr.rate.g2:

  Departure or exit rate for group 2 recovered persons. This parameter
  is only used for `SIR` model types.

- ...:

  Additional arguments passed to model.

## Value

An `EpiModel` object of class `param.net`.

## Details

`param.net` sets the epidemic parameters for the stochastic network
models simulated with the
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
function. Models may use the base types, for which these parameters are
used, or new process modules which may use these parameters (but not
necessarily). A detailed description of network model parameterization
for base models is found in the [Network Modeling for
Epidemics](https://epimodel.github.io/sismid/) tutorials.

For base models, the model specification will be chosen as a result of
the model parameters entered here and the control settings in
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md).
One-group and two-group models are available, where the latter assumes a
heterogeneous mixing between two distinct partitions in the population
(e.g., men and women). Specifying any two-group parameters (those with a
`.g2`) implies the simulation of a two-group model. All the parameters
for a desired model type must be specified, even if they are zero.

## The `act.rate` Parameter

A key difference between these network models and DCM/ICM classes is the
treatment of transmission events. With DCM and ICM, contacts or
partnerships are mathematically instantaneous events: they have no
duration in time, and thus no changes may occur within them over time.
In contrast, network models allow for partnership durations defined by
the dynamic network model, summarized in the model dissolution
coefficients calculated in
[`dissolution_coefs()`](https://epimodel.github.io/EpiModel/reference/dissolution_coefs.md).
Therefore, the `act.rate` parameter has a different interpretation here,
where it is the number of transmissible acts *per partnership* per unit
time.

## Parameters by Duration of Infection

The `inf.prob`, `act.rate`, and `rec.rate` arguments (and their `.g2`
companions) may vary with the duration of infection, given as
[`by_infection_duration()`](https://epimodel.github.io/EpiModel/reference/by_infection_duration.md)
objects, by stage or one value per time step since infection. For
example,
`inf.prob = by_infection_duration(c(acute = 0.5, chronic = 0.1), durations = c(2, Inf))`
gives a 0.5 transmission probability for the first two time steps of the
infected partner's infection and 0.1 from the third time step on, until
the person recovers, departs, or the simulation ends. This is variation
over the course of each infection, not over calendar time; to change a
parameter at a given time step of the simulation, use scenarios or
parameter updaters (see
[`vignette("model-parameters", package = "EpiModel")`](https://epimodel.github.io/EpiModel/articles/model-parameters.md)).

With the built-in model types, a plain vector of length greater than one
passed to these arguments is an error. Before EpiModel 2.7.0, the
built-in modules read such a vector as varying with the duration of
infection; wrap it in
[`by_infection_duration()`](https://epimodel.github.io/EpiModel/reference/by_infection_duration.md)
for the same model. Custom modules read their parameters as they choose,
so plain vectors remain available to them.

## Multi-Layer Parameters

In models with more than one network layer (a list of
[`netest()`](https://epimodel.github.io/EpiModel/reference/netest.md)
and
[`netclique()`](https://epimodel.github.io/EpiModel/reference/netclique.md)
objects passed to
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)),
the `inf.prob`, `inf.prob.g2`, and `act.rate` arguments may be given per
layer as
[`multilayer()`](https://epimodel.github.io/EpiModel/reference/multilayer.md)
objects with one entry per layer, in the order of the layer list:
`inf.prob = multilayer(0.45, 0.10)` sets a per-act transmission
probability of 0.45 on the first layer and 0.10 on the second. Each
entry may itself be a
[`by_infection_duration()`](https://epimodel.github.io/EpiModel/reference/by_infection_duration.md)
object. The built-in infection modules look up the entry for the layer
each discordant edge belongs to; a parameter that is not a `multilayer`
object applies to every layer. The transmission matrix records the layer
of each transmission in its `network` column.

## Using a Parameter data.frame

It is possible to set input parameters using a specifically formatted
`data.frame` object. The first 3 columns of this `data.frame` must be:

- `param`: The name of the parameter. If this is a non-scalar parameter
  (a vector of length \> 1), end the parameter name with the position on
  the vector (e.g., `"p_1"`, `"p_2"`, ...).

- `value`: the value for the parameter (or the value of the parameter in
  the Nth position if non-scalar).

- `type`: a character string containing either `"numeric"`, `"logical"`,
  or `"character"` to define the parameter object class.

In addition to these 3 columns, the `data.frame` can contain any number
of other columns, such as `details` or `source` columns to document
parameter meta-data. However, these extra columns will not be used by
EpiModel.

This data.frame is then passed in to `param.net` under a
`data.frame.params` argument. Further details and examples are provided
in the "Working with Model Parameters in EpiModel" vignette. Note that
releases through v2.6.1 documented this argument as
`data.frame.parameters`, which was never an accepted name; passing it
now produces an error pointing to `data.frame.params`.

## Parameter Uncertainty

To vary parameter values across simulations for uncertainty or
sensitivity analysis, draw the values in advance into a `data.frame`
with one row per draw, convert it with
[`create_scenario_list()`](https://epimodel.github.io/EpiModel/reference/create_scenario_list.md),
and apply each resulting scenario to the base parameters with
[`use_scenario()`](https://epimodel.github.io/EpiModel/reference/use_scenario.md).
The "Working with Model Parameters in EpiModel" vignette works through
an example. The `random.params` argument that served this purpose
through v2.6.2 has been removed; passing it now produces an error.

## Parameters with New Modules

To build original models outside of the base models, new process modules
may be constructed to replace the existing modules or to supplement the
existing set. These are passed into the control settings in
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md).
New modules may use either the existing model parameters named here, an
original set of parameters, or a combination of both. The `...` allows
the user to pass an arbitrary set of original model parameters into
`param.net`. Whereas there are strict checks with default modules for
parameter validity, this becomes a user responsibility when using new
modules.

## See also

Use
[`init.net()`](https://epimodel.github.io/EpiModel/reference/init.net.md)
to specify the initial conditions and
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md)
to specify the control settings. Run the parameterized model with
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md).

## Examples

``` r
# \donttest{
## Example SIR model parameterization
# Network model estimation
nw <- network_initialize(n = 100)
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#> Starting simulated annealing (SAN)
#> Iteration 1 of at most 4
#> Finished simulated annealing
#> Starting maximum pseudolikelihood estimation (MPLE):
#> Obtaining the responsible dyads.
#> Evaluating the predictor and response matrix.
#> Maximizing the pseudolikelihood.
#> Finished MPLE.

# Parameters, initial conditions, and control settings
param <- param.net(inf.prob = 0.3, act.rate = 2, rec.rate = 0.02)
param
#> Model Parameters
#> ---------------------------
#> inf.prob = 0.3
#> act.rate = 2
#> rec.rate = 0.02

init <- init.net(i.num = 10, r.num = 0)
control <- control.net(type = "SIR", nsteps = 10, nsims = 3, verbose = FALSE)

# Simulate the model
sim <- netsim(est, param, init, control)
sim
#> EpiModel Simulation
#> =======================
#> Model class: netsim
#> 
#> Simulation Summary
#> -----------------------
#> Model type: SIR
#> No. simulations: 3
#> No. time steps: 10
#> No. NW groups: 1
#> 
#> Model Parameters
#> ---------------------------
#> inf.prob = 0.3
#> act.rate = 2
#> rec.rate = 0.02
#> groups = 1
#> 
#> Model Output
#> -----------------------
#> Variables: s.num i.num r.num num si.flow ir.flow
#> Networks: sim1 ... sim3
#> Transmissions: sim1 ... sim3
#> 
#> Formation Statistics
#> ----------------------- 
#>       Target Sim Mean Pct Diff Sim SE Z Score SD(Sim Means) SD(Statistic)
#> edges     50       48       -4  0.826   -2.42         4.927         4.526
#> 
#> 
#> Duration Statistics
#> ----------------------- 
#>       Target Sim Mean Pct Diff Sim SE Z Score SD(Sim Means) SD(Statistic)
#> edges     20   19.299   -3.504  0.181  -3.867         0.426         0.991
#> 
#> Dissolution Statistics
#> ----------------------- 
#>       Target Sim Mean Pct Diff Sim SE Z Score SD(Sim Means) SD(Statistic)
#> edges   0.05    0.055     9.38  0.005   1.022         0.012         0.032
#> 
# }
```
