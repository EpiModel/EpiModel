# Parameters That Vary With the Duration of Infection

Marks transmission probabilities, act rates, or recovery rates as
varying with the time since infection, for the `inf.prob`,
`inf.prob.g2`, `act.rate`, `rec.rate`, and `rec.rate.g2` arguments of
[`param.net()`](https://epimodel.github.io/EpiModel/reference/param.net.md).
The values are given either by stage, with the length of each stage in
time steps, or one per time step since infection. The built-in
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
modules apply the value at each infected node's current duration of
infection.

## Usage

``` r
by_infection_duration(x, durations = NULL)
```

## Arguments

- x:

  A numeric vector of values. With `durations`, one value per stage of
  infection, optionally named for the stages (such as
  `c(acute = 0.25, chronic = 0.02)`). Without `durations`, one value per
  time step since infection: the first element applies in the first time
  step of an infection, the second in the second, and so on, and the
  last element carries forward for the rest of the infection.

- durations:

  Optional numeric vector of the same length as `x`, giving the length
  of each stage in time steps. Every stage but the last must last a
  whole number of time steps, one or more. The last stage lasts until
  the infection ends, so its duration must be `Inf`; for a value that
  changes after a finite last stage, add a stage with the value that
  follows.

## Value

An object of class `by_infection_duration`: the numeric vector of values
`x`, with the stage durations in its `durations` attribute. Without
`durations`, every stage lasts one time step but the last.

## Details

A duration of infection is the number of time steps since the infected
node was infected, counted from 1: a node infected at time step 10 is in
its first step of infection at time steps 10 and 11, its second at time
step 12, and so on. For transmission, the duration is that of the
infected partner; for recovery, that of the recovering node.

The two forms describe the same thing.
`by_infection_duration(c(acute = 0.25, chronic = 0.02), durations = c(10, Inf))`
gives 0.25 in the first ten time steps of an infection and 0.02 from the
eleventh on, which can represent an acute stage of high infectiousness;
`by_infection_duration(c(rep(0.25, 10), 0.02))` is the same parameter
given per time step. The stage form suits parameters reported by stage
of infection; the per-step form suits a profile computed elsewhere, such
as an infectiousness curve.

This is variation over the course of each infection, not over calendar
time: to change a parameter at a given time step of the simulation, use
scenarios or parameter updaters (see
[`vignette("model-parameters", package = "EpiModel")`](https://epimodel.github.io/EpiModel/articles/model-parameters.md)).

In multi-layer models, an entry of a
[`multilayer()`](https://epimodel.github.io/EpiModel/reference/multilayer.md)
parameter may itself be a `by_infection_duration` object:
`inf.prob = multilayer(by_infection_duration(c(0.5, 0.1)), 0.05)`.
Arithmetic on the object, such as `x * 2`, scales the values and keeps
the stages.

With the built-in modules, a plain vector passed to these arguments
stops
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
with an error, since a plain vector carries no record of what its
positions mean; custom modules keep their own reading of plain vectors,
such as values by group or by layer. Parameter tables cannot represent
these objects yet, so
[`param.net_to_table()`](https://epimodel.github.io/EpiModel/reference/param.net_to_table.md)
stops on them.

## See also

[`param.net()`](https://epimodel.github.io/EpiModel/reference/param.net.md),
[`multilayer()`](https://epimodel.github.io/EpiModel/reference/multilayer.md).

## Examples

``` r
# An acute stage of ten time steps with a higher transmission probability
inf.prob <- by_infection_duration(c(acute = 0.25, chronic = 0.02),
                                  durations = c(10, Inf))
inf.prob
#> Values by stage of infection (time steps since infection):
#>    Stage Steps Value
#>    acute  1-10  0.25
#>  chronic   11+  0.02

# Recovery impossible in the first 20 time steps of an infection, then
# certain
by_infection_duration(c(0, 1), durations = c(20, Inf))
#> Values by stage of infection (time steps since infection):
#>  Stage Steps Value
#>      1  1-20     0
#>      2   21+     1

# One value per time step, here an infectiousness profile
by_infection_duration(round(dgamma(1:15, shape = 3, rate = 0.5), 3))
#> Values by duration of infection (time steps since infection), the last carried forward:
#>  [1] 0.038 0.092 0.126 0.135 0.128 0.112 0.092 0.073 0.056 0.042 0.031 0.022
#> [13] 0.016 0.011 0.008

param <- param.net(inf.prob = inf.prob, act.rate = 1)
param
#> Model Parameters
#> ---------------------------
#> inf.prob = 
#> by_infection_duration(c(acute = 0.25, chronic = 0.02), durations = c(10, Inf))
#> act.rate = 1

if (FALSE) { # \dontrun{
nw <- network_initialize(n = 100)
est <- netest(nw, formation = ~edges, target.stats = 50,
              coef.diss = dissolution_coefs(~offset(edges), 10),
              verbose = FALSE)
sim <- netsim(est, param, init.net(i.num = 10),
              control.net(type = "SI", nsteps = 25, nsims = 1,
                          verbose = FALSE))
tm <- get_transmat(sim)
table(tm$infDur, tm$transProb)
} # }
```
