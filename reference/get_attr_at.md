# Get Tracked Attributes at a Given Time Step

Reconstructs the state of all tracked attributes at a specific time step
from the delta-compressed attribute history. Returns a `tibble` with one
row per active node and one column per tracked attribute.

## Usage

``` r
get_attr_at(sim, at, sim_num = 1, compute_age = FALSE)
```

## Arguments

- sim:

  An `EpiModel` object of class `netsim`.

- at:

  The time step at which to reconstruct attributes, on the time step
  labels of `sim` (see
  [`truncate_sim()`](https://epimodel.github.io/EpiModel/reference/truncate_sim.md)).
  It must lie between the start of the attribute history and the last
  time step. The history starts at the first time step for a new
  simulation, and at `start - 1` for a restarted one, whose history
  starts again from a full snapshot of the state it restarted from.

- sim_num:

  The simulation number to use (default 1).

- compute_age:

  If `TRUE` and `age` is tracked, compute the age of nodes at `at` from
  their last recorded age, assuming age in years and 364 days in a year.
  This needs a `time.unit` parameter (the number of days in a time step)
  in `sim$param`, passed to
  [`param.net`](https://epimodel.github.io/EpiModel/reference/param.net.md)
  or set by the initialization module (as EpiModelHIV does from its
  `epistats`).

## Value

A `tibble` with columns `unique_id`, `active`, and one column per
tracked attribute, containing the most recent value at or before `at`
for each active node.

## See also

[`record_attr_history()`](https://epimodel.github.io/EpiModel/reference/record_attr_history.md),
[`get_attr_history()`](https://epimodel.github.io/EpiModel/reference/get_attr_history.md)

## Examples

``` r
# \donttest{
nw <- network_initialize(n = 100)
est <- netest(nw, formation = ~edges, target.stats = 50,
              coef.diss = dissolution_coefs(~offset(edges), 10),
              verbose = FALSE)
#> Starting simulated annealing (SAN)
#> Iteration 1 of at most 4
#> Finished simulated annealing
#> Starting maximum pseudolikelihood estimation (MPLE):
#> Obtaining the responsible dyads.
#> Evaluating the predictor and response matrix.
#> Maximizing the pseudolikelihood.
#> Finished MPLE.
control <- control.net(type = "SI", nsteps = 20, nsims = 1,
                       tracked.attributes = "status", verbose = FALSE)
sim <- netsim(est, param.net(inf.prob = 0.3, act.rate = 1),
              init.net(i.num = 10), control)

# The status of every node at time step 10
at10 <- get_attr_at(sim, at = 10)
head(at10)
#> # A tibble: 6 × 3
#>   unique_id active status
#>       <int>  <int> <chr> 
#> 1         1      1 i     
#> 2         2      1 i     
#> 3         3      1 i     
#> 4         4      1 s     
#> 5         5      1 s     
#> 6         6      1 s     
table(at10$status)
#> 
#>  i  s 
#> 36 64 
# }
```
