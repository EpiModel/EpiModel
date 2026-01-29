# Summary Model Statistics

Extracts and prints model statistics simulated with `netsim`.

## Usage

``` r
# S3 method for class 'netsim'
summary(object, at, digits = 3, ...)
```

## Arguments

- object:

  An `EpiModel` object of class `netsim`.

- at:

  Time step for model statistics.

- digits:

  Number of significant digits to print.

- ...:

  Additional summary function arguments.

## Details

This function provides summary statistics for the main epidemiological
outcomes (state and transition size and prevalence) from a `netsim`
model. Time-specific summary measures are provided, so it is necessary
to input a time of interest.

## See also

[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## SI Model without Network Feedback
# Initialize network and set network model parameters
nw <- network_initialize(n = 100)
nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)

# Estimate the ERGM models (see help for netest)
est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

# Parameters, initial conditions, and controls for model
param <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.15)
init <- init.net(i.num = 10, i.num.g2 = 10)
control <- control.net(type = "SI", nsteps = 100, nsims = 5, verbose.int = 0)

# Run the model simulation
mod <- netsim(est1, param, init, control)

summary(mod, at = 1)
summary(mod, at = 50)
summary(mod, at = 100)
} # }
```
