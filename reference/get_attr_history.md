# Extract the Attributes History from Network Simulations

Extract the Attributes History from Network Simulations

## Usage

``` r
get_attr_history(sims)
```

## Arguments

- sims:

  An `EpiModel` object of class `netsim`.

## Value

A list of `data.frame`s, one for each "measure" recorded in the
simulation by the `record_attr_history` function.

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

# One data frame per tracked attribute, with the changes at each step
hist <- get_attr_history(sim)
names(hist)
#> [1] "active" "status"
head(hist$status)
#> # A tibble: 6 × 5
#>     sim  time attribute  uids values
#>   <dbl> <int> <chr>     <int> <chr> 
#> 1     1     1 status        1 i     
#> 2     1     1 status        2 s     
#> 3     1     1 status        3 s     
#> 4     1     1 status        4 s     
#> 5     1     1 status        5 s     
#> 6     1     1 status        6 s     
# }
```
