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
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

- sim:

  Simulation number of extracted network.

- deduplicate:

  If `TRUE`, randomly select one transmission event in the case that
  multiple events current per newly infected agent within a time step.

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
## Simulate SI epidemic on two-group Bernoulli random graph
nw <- network_initialize(n = 100)
nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
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
#> # A tibble: 20 × 8
#> # Groups:   at, sus [20]
#>       at   sus   inf network infDur transProb actRate finalProb
#>    <int> <int> <int>   <int>  <dbl>     <dbl>   <dbl>     <dbl>
#>  1     2     1    43       1      7      0.3        1      0.3 
#>  2     2    12    13       1      6      0.3        1      0.3 
#>  3     2    19    34       1      1      0.3        1      0.3 
#>  4     3     2    57       1      7      0.3        1      0.3 
#>  5     3    47    13       1      7      0.3        1      0.3 
#>  6     3    64    12       1      1      0.15       1      0.15
#>  7     3    99    97       1      5      0.15       1      0.15
#>  8     4    25     2       1      1      0.3        1      0.3 
#>  9     4    62    64       1      1      0.15       1      0.15
#> 10     4    73    58       1      3      0.15       1      0.15
#> 11     5    63    25       1      1      0.15       1      0.15
#> 12     5    68    67       1      5      0.15       1      0.15
#> 13     6    82    31       1      6      0.15       1      0.15
#> 14     7    60    68       1      2      0.15       1      0.15
#> 15     7    61    68       1      2      0.15       1      0.15
#> 16     8     9     1       1      6      0.3        1      0.3 
#> 17     8    86    19       1      6      0.15       1      0.15
#> 18     9    27    64       1      6      0.3        1      0.3 
#> 19     9    74    31       1      9      0.15       1      0.15
#> 20    10    83    31       1     10      0.15       1      0.15
```
