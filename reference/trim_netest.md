# Function to Reduce the Size of a `netest` Object

Trims formula environments from the `netest` object. Optionally converts
the `newnetwork` element of the `netest` object to a `networkLite`
class, and removes the `fit` element (if present) from the `netest`
object.

## Usage

``` r
trim_netest(
  object,
  as.networkLite = TRUE,
  keep.fit = FALSE,
  keep = character(0)
)
```

## Arguments

- object:

  A `netest` class object.

- as.networkLite:

  If `TRUE`, converts `object$newnetwork` to a `networkLite`.

- keep.fit:

  If `FALSE`, removes the `object$fit` (if present) on the `netest`
  object.

- keep:

  Character vector of object names to keep in formula environments. By
  default, all objects are removed.

## Value

A `netest` object with formula environments trimmed, optionally with the
`newnetwork` element converted to a `networkLite` and the `fit` element
removed.

## Details

With larger, more complex network structures with epidemic models, it is
generally useful to reduce the memory footprint of the fitted TERGM
model object (estimated with
[`netest()`](http://epimodel.github.io/EpiModel/reference/netest.md)).
This utility function removes all but the bare essentials needed for
simulating a network model with
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

The function always trims the environments of `object$constraints` and
`object$coef.diss$dissolution`.

When both `edapprox = TRUE` and `nested.edapprox = TRUE` in the `netest`
call, also trims the environments of `object$formula` and
`object$formation`.

When both `edapprox = TRUE` and `nested.edapprox = FALSE` in the
`netest` call, also trims the environments of `object$formula`,
`environment(object$formation)$formation`, and
`environment(object$formation)$dissolution`.

When `edapprox = FALSE` in the `netest` call, also trims the
environments of `object$formation`,
`environment(object$formula)$formation` and
`environment(object$formula)$dissolution`.

By default all objects are removed from these trimmed environments.
Specific objects may be retained by passing their names as the `keep`
argument. For the output of `trim_netest` to be usable in
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md)
simulation, any objects referenced in the formulas should be included in
the `keep` argument.

If `as.networkLite = TRUE`, converts `object$newnetwork` to a
`networkLite` object. If `keep.fit = FALSE`, removes `fit` (if present)
from `object`.

## Examples

``` r
nw <- network_initialize(n = 100)
formation <- ~edges + concurrent
target.stats <- c(50, 25)
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
est <- netest(nw, formation, target.stats, coef.diss,
              set.control.ergm = control.ergm(MCMC.burnin = 1e5,
                                              MCMC.interval = 1000))
#> Starting maximum pseudolikelihood estimation (MPLE):
#> Obtaining the responsible dyads.
#> Evaluating the predictor and response matrix.
#> Maximizing the pseudolikelihood.
#> Finished MPLE.
#> Starting Monte Carlo maximum likelihood estimation (MCMLE):
#> Iteration 1 of at most 60:
#> 1 
#> Optimizing with step length 1.0000.
#> The log-likelihood improved by 0.3250.
#> Estimating equations are not within tolerance region.
#> Iteration 2 of at most 60:
#> 1 
#> Optimizing with step length 1.0000.
#> The log-likelihood improved by 0.0237.
#> Convergence test p-value: < 0.0001. 
#> Converged with 99% confidence.
#> Finished MCMLE.
#> This model was fit using MCMC.  To examine model diagnostics and check
#> for degeneracy, use the mcmc.diagnostics() function.
print(object.size(est), units = "KB")
#> 423.5 Kb

est.small <- trim_netest(est)
print(object.size(est.small), units = "KB")
#> 58.8 Kb
```
