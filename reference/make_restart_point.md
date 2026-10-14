# Make a Lightweight Restart Point From a `netsim` Object with tergmLite

Extract the elements required for re-initializing a `netsim` simulation
from a completed simulation. This function also resets the Unique IDs
and Time values to reduce the size of the simulation. This function only
works for simulations run with `tergmLite = TRUE` and `save.run = TRUE`
in
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md).

## Usage

``` r
make_restart_point(sim_obj, time_attrs, sim_num = 1, keep_steps = 1)
```

## Arguments

- sim_obj:

  a `netsim` object from an ended `netsim` call, run with
  `tergmLite = TRUE` and `save.run = TRUE`.

- time_attrs:

  a `character` vector containing the names of the attributes that are
  expressed in time-steps. These will be offset so the last step in the
  original simulation become the step 1 (default) in the new ones. The
  core attributes `entrTime`, `exitTime`, and `infTime` (when present)
  are always offset. If the model has no other such attributes, pass
  [`c()`](https://rdrr.io/r/base/c.html).

- sim_num:

  the number of the simulation to extract from the `netsim` object
  (default = 1).

- keep_steps:

  The number of simulation steps to keep from the previous run. By
  default only keep one but more is possible if some back-history is
  wanted.

## Value

a trimmed `netsim` object with only one simulation that is ready to be
used as a restart point.

## Details

The restart point created always contains a single simulation and drops
the `attr.history` and the `raw.records` from the initial simulation.

`sim_obj` does not need to hold the parameters: when restarting,
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
takes them from its `param` argument alone, which must hold all the
parameters of the model.

The epi trackers, cumulative edgelists, transmission matrix and
`nwstats` are truncated to only contain the last `keep_steps` entries.
The unique IDs in the transmission matrix are re-based like the ones of
the attributes.

Warning: the `time_attrs` argument is mandatory. Almost all extension
models have such attributes (e.g. the time of HIV diagnosis). If no such
attribute exists beyond the core ones, passing
[`c()`](https://rdrr.io/r/base/c.html) will allow the function to run
while ensuring that this was done on purpose. An attribute holding a
time step that is not offset puts the events it records in the future of
the restarted simulation: for example, a duration of infection computed
from an `infTime` that was not offset would be negative.

A dynamic
[`netcensus()`](https://epimodel.github.io/EpiModel/reference/netcensus.md)
layer reads the observed network on the simulation clock, which this
function re-bases, so it stops on such a layer. Restart those
simulations from the `netsim` object itself, selecting a simulation with
[`get_sims()`](https://epimodel.github.io/EpiModel/reference/get_sims.md).

When restarting from the output of this function, it is suggested to
express the time steps in a relative manner in `control.net`:

    control.net(
      start = restart_point$control$nsteps + 1,
      nsteps = restart_point$control$nsteps + 1 + 104
    )

## Examples

``` r
# \donttest{
nw <- network_initialize(n = 100)
est <- netest(nw, formation = ~edges, target.stats = 50,
              coef.diss = dissolution_coefs(~offset(edges), 20),
              verbose = FALSE)
#> Starting simulated annealing (SAN)
#> Iteration 1 of at most 4
#> Finished simulated annealing
#> Starting maximum pseudolikelihood estimation (MPLE):
#> Obtaining the responsible dyads.
#> Evaluating the predictor and response matrix.
#> Maximizing the pseudolikelihood.
#> Finished MPLE.
param <- param.net(inf.prob = 0.3, act.rate = 1)
init <- init.net(i.num = 10)
control <- control.net(type = "SI", nsteps = 20, nsims = 2,
                       tergmLite = TRUE, resimulate.network = TRUE,
                       save.run = TRUE, verbose = FALSE)
sim <- netsim(est, param, init, control)

# Make a restart point from simulation 2. The built-in SI model has no time
# attribute beyond the core ones, which are offset by default.
x <- make_restart_point(sim, time_attrs = c(), sim_num = 2)

# Run it for 10 more time steps
control <- control.net(type = "SI", start = x$control$nsteps + 1,
                       nsteps = x$control$nsteps + 10,
                       tergmLite = TRUE, resimulate.network = TRUE,
                       verbose = FALSE)
sim2 <- netsim(x, param, init, control)
# }
```
