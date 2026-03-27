# Extract Transmissions Matrix from Network Epidemic Model

Extracts the matrix of transmission data for each transmission event
that occurred within a network epidemic model.

## Usage

``` r
is.transmat(x)

get_transmat(x, sim = 1, deduplicate = TRUE)
```

## Arguments

- x:

  An `EpiModel` object of class
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md).

- sim:

  Simulation number of extracted network.

- deduplicate:

  If `TRUE`, randomly select one transmission event in the case that
  multiple events occur per newly infected agent within a time step.

## Value

A data frame with the following standard columns:

- **at:** the time step at which the transmission occurred.

- **sus:** the ID number of the susceptible (newly infected) node.

- **inf:** the ID number of the infecting node.

- **infDur:** the duration of the infecting node's disease at the time
  of the transmission.

- **transProb:** the probability of transmission per act.

- **actRate:** the rate of acts per unit time.

- **finalProb:** the final transmission probability for the transmission
  event.

## Examples

``` r
# \donttest{
## Simulate SI epidemic on two-group Bernoulli random graph
nw <- network_initialize(n = 100)
nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#> Starting simulated annealing (SAN)
#> Iteration 1 of at most 4
#> Finished simulated annealing
#> Starting maximum pseudolikelihood estimation (MPLE):
#> Obtaining the responsible dyads.
#> Evaluating the predictor and response matrix.
#> Maximizing the pseudolikelihood.
#> Finished MPLE.
param <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.15)
init <- init.net(i.num = 10, i.num.g2 = 10)
control <- control.net(type = "SI", nsteps = 10, nsims = 3, verbose = FALSE)
mod <- netsim(est, param, init, control)

## Extract the transmission matrix from simulation 2
get_transmat(mod, sim = 2)
#> # A tibble: 24 × 8
#> # Groups:   at, sus [24]
#>       at   sus   inf network infDur transProb actRate finalProb
#>    <int> <int> <int>   <int>  <dbl>     <dbl>   <dbl>     <dbl>
#>  1     2     1    68       1      9      0.3        1      0.3 
#>  2     2    19    81       1      8      0.3        1      0.3 
#>  3     2    21    92       1      8      0.3        1      0.3 
#>  4     2    43    77       1      9      0.3        1      0.3 
#>  5     3    13    81       1      9      0.3        1      0.3 
#>  6     3    35     8       1      9      0.3        1      0.3 
#>  7     3    36    42       1      3      0.3        1      0.3 
#>  8     4    20    16       1      7      0.3        1      0.3 
#>  9     4    45    78       1      5      0.3        1      0.3 
#> 10     4    84    36       1      1      0.15       1      0.15
#> # ℹ 14 more rows
# }
```
