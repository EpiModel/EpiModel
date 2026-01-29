# Extract Network Simulations

Subsets the entire `netsim` object to a subset of simulations,
essentially functioning like a reverse of `merge`.

## Usage

``` r
get_sims(x, sims = NULL, var = NULL)
```

## Arguments

- x:

  An object of class `netsim`.

- sims:

  Either a numeric vector of simulation numbers to retain in the output
  object, or `"mean"`, which selects the one simulation with the value
  of the variable specified in `var` closest to the mean of `var` across
  all simulations.

- var:

  A character vector of variables to retain from `x` if `sims` is a
  numeric vector, or a single variable name for selecting the average
  simulation from the set if `sims = "mean"`.

## Value

An updated object of class `netsim` containing only the simulations
specified in `sims` and the variables specified in `var`.

## Examples

``` r
# Network model estimation
nw <- network_initialize(n = 100)
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#> Starting maximum pseudolikelihood estimation (MPLE):
#> Obtaining the responsible dyads.
#> Evaluating the predictor and response matrix.
#> Maximizing the pseudolikelihood.
#> Finished MPLE.

# Epidemic model
param <- param.net(inf.prob = 0.3)
init <- init.net(i.num = 10)
control <- control.net(type = "SI", nsteps = 10, nsims = 3, verbose.int = 0)
mod1 <- netsim(est1, param, init, control)
#> 
#> Starting Network Simulation...
#> Sim = 1/3
#> Sim = 2/3
#> Sim = 3/3

# Get sim 2
s.g2 <- get_sims(mod1, sims = 2)

# Get sims 2 and 3 and keep only a subset of variables
s.g2.small <- get_sims(mod1, sims = 2:3, var = c("i.num", "si.flow"))

# Extract the mean simulation for the variable i.num
sim.mean <- get_sims(mod1, sims = "mean", var = "i.num")
```
