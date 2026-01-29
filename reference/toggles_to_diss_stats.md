# Convert Matrix of Toggles to Dissolution and Duration Statistics

Convert Matrix of Toggles to Dissolution and Duration Statistics

## Usage

``` r
toggles_to_diss_stats(toggles, coef.diss, nsteps, nw, time.start = 0L)
```

## Arguments

- toggles:

  A matrix of toggles, as produced by
  [`tedgelist_to_toggles()`](http://epimodel.github.io/EpiModel/reference/tedgelist_to_toggles.md).

- coef.diss:

  Dissolution coefficients used in the simulation.

- nsteps:

  Number of time steps in the simulation.

- nw:

  Network used in the simulation.

- time.start:

  Starting time for the simulation.

## Value

Named list containing dissolution and duration statistics matrices and
other related information.
