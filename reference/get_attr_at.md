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

  The time step at which to reconstruct attributes.

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
