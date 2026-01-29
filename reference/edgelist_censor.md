# Table of Edge Censoring

Outputs a table of the number and percent of edges that are
left-censored, right-censored, both-censored, or uncensored for a
`networkDynamic` object.

## Usage

``` r
edgelist_censor(el)
```

## Arguments

- el:

  A timed edgelist with start and end times extracted from a
  `networkDynamic` object using the `as.data.frame.networkDynamic`
  function.

## Value

A 4 x 2 table containing the number and percent of edges in `el` that
are left-censored, right-censored, both-censored, or uncensored.

## Details

Given a STERGM simulation over a specified number of time steps, the
edges within that simulation may be left-censored (started before the
first step), right-censored (continued after the last step), right and
left-censored, or uncensored. The amount of censoring will increase when
the average edge duration approaches the length of the simulation.

## Examples

``` r
# Initialize and parameterize network model
nw <- network_initialize(n = 100)
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)

# Model estimation
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#> Starting maximum pseudolikelihood estimation (MPLE):
#> Obtaining the responsible dyads.
#> Evaluating the predictor and response matrix.
#> Maximizing the pseudolikelihood.
#> Finished MPLE.

# Simulate the network and extract a timed edgelist
dx <- netdx(est, nsims = 1, nsteps = 100, keep.tedgelist = TRUE,
      verbose = FALSE)
#> Warning: NAs introduced by coercion to integer range
el <- as.data.frame(dx)

# Calculate censoring
edgelist_censor(el)
#>             num      pct
#> Left Cens.    0 0.000000
#> Right Cens.  53 0.171521
#> Both Cens.    0 0.000000
#> No Cens.    256 0.828479
```
