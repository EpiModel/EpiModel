# Create a Summary Table of Simulation Statistics

Create a Summary Table of Simulation Statistics

## Usage

``` r
make_stats_table(stats, targets)
```

## Arguments

- stats:

  A list of simulated statistics matrices, of length equal to the number
  of simulations performed. Each matrix should have one row for each
  simulated network if `dynamic == FALSE`, one row for each time step if
  `dynamic == TRUE`, and one column for each statistic. The columns
  should be named for the statistics they correspond to, with all
  matrices having the same statistics, in the same order. Each matrix
  may have an `attr`-style attribute named `"ess"` attached, giving the
  effective sample sizes for the columns of the matrix; if this
  attribute is `NULL`, then the effective sample sizes will be computed
  within the call to `make_stats_table`.

- targets:

  A vector of target values for the statistics in `stats`. May be named
  (in which case targets will be matched to statistics based on column
  names in matrices in `stats`) or unnamed (in which case targets will
  be matched to statistics based on position, and the number of targets
  must equal the number of columns).

## Value

A `data.frame` summarizing the simulated statistics.
