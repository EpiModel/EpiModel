# Epidemic Parameters for Stochastic Network Models

Sets the epidemic parameters for stochastic network models simulated
with
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

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
  of infection to the group 1 nodes. This may also be a vector of
  probabilities, with each element corresponding to the probability in
  that time step of infection (see Time-Varying Parameters below).

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
  (see `act.rate` Parameter below). This may also be a vector of rates,
  with each element corresponding to the rate in that time step of
  infection (see Time-Varying Parameters below).

- rec.rate:

  Average rate of recovery with immunity (in `SIR` models) or
  re-susceptibility (in `SIS` models). The recovery rate is the
  reciprocal of the disease duration. For two-group models, this is the
  recovery rate for group 1 persons only. This parameter is only used
  for `SIR` and `SIS` models. This may also be a vector of rates, with
  each element corresponding to the rate in that time step of infection
  (see Time-Varying Parameters below).

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
  probability of transmission to group 2 members.

- rec.rate.g2:

  Average rate of recovery with immunity (in `SIR` models) or
  re-susceptibility (in `SIS` models) for group 2 persons. This
  parameter is only used for two-group `SIR` and `SIS` models.

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
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md)
function. Models may use the base types, for which these parameters are
used, or new process modules which may use these parameters (but not
necessarily). A detailed description of network model parameterization
for base models is found in the [Network Modeling for
Epidemics](https://epimodel.github.io/sismid/) tutorials.

For base models, the model specification will be chosen as a result of
the model parameters entered here and the control settings in
[`control.net()`](http://epimodel.github.io/EpiModel/reference/control.net.md).
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
[`dissolution_coefs()`](http://epimodel.github.io/EpiModel/reference/dissolution_coefs.md).
Therefore, the `act.rate` parameter has a different interpretation here,
where it is the number of transmissible acts *per partnership* per unit
time.

## Time-Varying Parameters

The `inf.prob`, `act.rate`, `rec.rate` arguments (and their `.g2`
companions) may be specified as time-varying parameters by passing in a
vector of probabilities or rates, respectively. The value in each
position on the vector then corresponds to the probability or rate at
that discrete time step for the infected partner. For example, an
`inf.prob` of `c(0.5, 0.5, 0.1)` would simulate a 0.5 transmission
probability for the first two time steps of a person's infection,
followed by a 0.1 for the third time step. If the infected person has
not recovered or exited the population by the fourth time step, the
third element in the vector will carry forward until one of those events
occurs or the simulation ends. For further examples, see the [Network
Modeling for Epidemics](https://epimodel.github.io/sismid/) tutorials.

## Random Parameters

In addition to deterministic parameters in either fixed or time-varying
varieties above, one may also include a generator for random parameters.
These might include a vector of potential parameter values or a
statistical distribution definition; in either case, one draw from the
generator would be completed per individual simulation. This is possible
by passing a list named `random.params` into `param.net`, with each
element of `random.params` a named generator function. See the help page
and examples in
[`generate_random_params()`](http://epimodel.github.io/EpiModel/reference/generate_random_params.md).
A simple factory function for sampling is provided with
[`param_random()`](http://epimodel.github.io/EpiModel/reference/param_random.md)
but any function will do.

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
`data.frame.parameters` argument. Further details and examples are
provided in the "Working with Model Parameters in EpiModel" vignette.

## Parameters with New Modules

To build original models outside of the base models, new process modules
may be constructed to replace the existing modules or to supplement the
existing set. These are passed into the control settings in
[`control.net()`](http://epimodel.github.io/EpiModel/reference/control.net.md).
New modules may use either the existing model parameters named here, an
original set of parameters, or a combination of both. The `...` allows
the user to pass an arbitrary set of original model parameters into
`param.net`. Whereas there are strict checks with default modules for
parameter validity, this becomes a user responsibility when using new
modules.

## See also

Use
[`init.net()`](http://epimodel.github.io/EpiModel/reference/init.net.md)
to specify the initial conditions and
[`control.net()`](http://epimodel.github.io/EpiModel/reference/control.net.md)
to specify the control settings. Run the parameterized model with
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

## Examples

``` r
## Example SIR model parameterization with fixed and random parameters
# Network model estimation
nw <- network_initialize(n = 100)
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#> Starting maximum pseudolikelihood estimation (MPLE):
#> Obtaining the responsible dyads.
#> Evaluating the predictor and response matrix.
#> Maximizing the pseudolikelihood.
#> Finished MPLE.

# Random epidemic parameter list (here act.rate values are sampled uniformly
# with helper function param_random, and inf.prob follows a general Beta
# distribution with the parameters shown below)
my_randoms <- list(
  act.rate = param_random(1:3),
  inf.prob = function() rbeta(1, 1, 2)
)

# Parameters, initial conditions, and control settings
param <- param.net(rec.rate = 0.02, random.params = my_randoms)

# Printing parameters shows both fixed and and random parameter functions
param
#> Fixed Parameters
#> ---------------------------
#> rec.rate = 0.02
#> 
#> Random Parameters
#> (Not drawn yet)
#> ---------------------------
#> act.rate = <function>
#> inf.prob = <function>

# Set initial conditions and controls
init <- init.net(i.num = 10, r.num = 0)
control <- control.net(type = "SIR", nsteps = 10, nsims = 3, verbose = FALSE)

# Simulate the model
sim <- netsim(est, param, init, control)

# Printing the sim object shows the randomly drawn values for each simulation
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
#> Fixed Parameters
#> ---------------------------
#> rec.rate = 0.02
#> groups = 1
#> 
#> Random Parameters
#> ---------------------------
#> act.rate = 2 1 1
#> inf.prob = 0.5602979 0.2499597 0.5157282
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
#> edges     50   51.333    2.667  1.314   1.015         7.088         6.354
#> 
#> 
#> Duration Statistics
#> ----------------------- 
#>       Target Sim Mean Pct Diff Sim SE Z Score SD(Sim Means) SD(Statistic)
#> edges     20   19.693   -1.535  0.371  -0.827         1.269         1.478
#> 
#> Dissolution Statistics
#> ----------------------- 
#>       Target Sim Mean Pct Diff Sim SE Z Score SD(Sim Means) SD(Statistic)
#> edges   0.05    0.054    8.722  0.005   0.942         0.002         0.025
#> 

# Parameter sets can be extracted with:
get_param_set(sim)
#>   sim rec.rate vital groups act.rate  inf.prob
#> 1   1     0.02 FALSE      1        2 0.5602979
#> 2   2     0.02 FALSE      1        1 0.2499597
#> 3   3     0.02 FALSE      1        1 0.5157282
```
