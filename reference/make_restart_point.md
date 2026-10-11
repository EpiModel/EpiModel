# Make a Lightweight Restart Point From a `netsim` Object with tergmLite

Extract the elements required for re-initializing a `netsim` simulation
from a completed simulation. This function also resets the Unique IDs
and Time values to reduce the size of the simulation. This function only
works for simulations where `control$tergmLite = TRUE`

## Usage

``` r
make_restart_point(sim_obj, time_attrs, sim_num = 1, keep_steps = 1)
```

## Arguments

- sim_obj:

  a `netsim` object from an ended `netsim` call.

- time_attrs:

  a `character` vector containing the names of the attributes that are
  expressed in time-steps. These will be offset so the last step in the
  original simulation become the step 1 (default) in the new ones. If no
  such attributes exist, pass [`c()`](https://rdrr.io/r/base/c.html).

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

Warning: the `time_attrs` argument is mandatory. Almost all simulations
worth restarting have such attributes (e.g. time.of.hiv.infection). If
no such argument exists, passing [`c()`](https://rdrr.io/r/base/c.html)
will allow the function to run while ensuring that this was done on
purpose.

When restarting from the output of this function, it is suggested to
express the time steps in a relative manner in `control.net`:

    control.net(
      start = restart_point$control$nsteps + 1,
      nsteps = restart_point$control$nsteps + 1 + 104
    )

## Examples

``` r
if (FALSE) { # \dontrun{
# With  pre-existing `sim`, `param` and `init` object (see `netsim`)

# List all attributes that store a time step
time_attrs <- c(
  "inf.time",
  "stage.time",
  "aids.time",
  "prep.start.last"
)
# Make a restart point from simulation 1, re-run for 10 more timesteps
x <- make_restart_point(sim, time_attrs, sim_num = 1, keep_steps = 1)
control <- control_msm(
  start = x$control$nsteps + 1,
  nsteps = x$control$nsteps + 1 + 10
)
sim <- netsim(x, param, init, control)
} # }
```
