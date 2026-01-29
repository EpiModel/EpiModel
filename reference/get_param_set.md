# Extract the Parameter Set from Network Simulations

Extract the Parameter Set from Network Simulations

## Usage

``` r
get_param_set(sims)
```

## Arguments

- sims:

  An `EpiModel` object of class `netsim`.

## Value

A `data.frame` with one row per simulation and one column per parameter
or parameter element where the parameters are of size \> 1.

## Output Format

The outputted `data.frame` has one row per simulation and the columns
correspond to the parameters used in this simulation.

The column name will match the parameter name if it is a size 1
parameter or if the parameter is of size \> 1, there will be N columns
(with N being the size of the parameter) named `parameter.name_1`,
`parameter.name_2`, ..., `parameter.name_N`.

## Examples

``` r
# Setup network
nw <- network_initialize(n = 50)

est <- netest(
  nw, formation = ~edges,
  target.stats = c(25),
  coef.diss = dissolution_coefs(~offset(edges), 10, 0),
  verbose = FALSE
)
#> Starting maximum pseudolikelihood estimation (MPLE):
#> Obtaining the responsible dyads.
#> Evaluating the predictor and response matrix.
#> Maximizing the pseudolikelihood.
#> Finished MPLE.

init <- init.net(i.num = 10)

n <- 5

related.param <- data.frame(
  dummy.param = rbeta(n, 1, 2)
)

 my.randoms <- list(
   act.rate = param_random(c(0.25, 0.5, 0.75)),
   dummy.param = function() rbeta(1, 1, 2),
   dummy.strat.param = function() c(
     rnorm(1, 0, 10),
     rnorm(1, 10, 1)
   )
 )

param <- param.net(
  inf.prob = 0.3,
  dummy = c(0, 1, 2),
  random.params = my.randoms
)

control <- control.net(type = "SI", nsims = 3, nsteps = 5, verbose = FALSE)
mod <- netsim(est, param, init, control)

get_param_set(mod)
#>   sim inf.prob dummy_1 dummy_2 dummy_3 vital groups act.rate dummy.param
#> 1   1      0.3       0       1       2 FALSE      1     0.50   0.1676224
#> 2   2      0.3       0       1       2 FALSE      1     0.50   0.2262098
#> 3   3      0.3       0       1       2 FALSE      1     0.75   0.1882213
#>   dummy.strat.param_1 dummy.strat.param_2
#> 1           12.960228            10.79399
#> 2           -2.332582            10.62762
#> 3           12.932489            10.81655
```
