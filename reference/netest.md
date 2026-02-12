# Dynamic Network Model Estimation

Estimates statistical network models using the exponential random graph
modeling (ERGM) framework with extensions for dynamic/temporal models
(STERGM).

## Usage

``` r
netest(
  nw,
  formation,
  target.stats,
  coef.diss,
  constraints = NULL,
  coef.form = NULL,
  edapprox = TRUE,
  set.control.ergm = control.ergm(),
  set.control.tergm = control.tergm(MCMC.maxchanges = Inf),
  set.control.ergm.ego = NULL,
  verbose = FALSE,
  nested.edapprox = TRUE,
  ...
)
```

## Arguments

- nw:

  An object of class `network` or `egor`, with the latter indicating an
  `ergm.ego` fit.

- formation:

  Right-hand sided STERGM formation formula in the form `~edges + ...`,
  where `...` are additional network statistics.

- target.stats:

  Vector of target statistics for the formation model, with one number
  for each network statistic in the model. Ignored if fitting via
  `ergm.ego`.

- coef.diss:

  An object of class `disscoef` output from the
  [`dissolution_coefs`](http://epimodel.github.io/EpiModel/reference/dissolution_coefs.md)
  function.

- constraints:

  Right-hand sided formula specifying constraints for the modeled
  network, in the form `~...`, where `...` are constraint terms. By
  default, no constraints are set.

- coef.form:

  Vector of coefficients for the offset terms in the formation formula.

- edapprox:

  If `TRUE`, use the indirect edges dissolution approximation method for
  the dynamic model fit, otherwise use the more time-intensive full
  STERGM estimation (see details). For `nw` of class `egor`, only
  `edapprox = TRUE` is supported.

- set.control.ergm:

  Control arguments passed to `ergm` (see details).

- set.control.tergm:

  Control arguments passed to `tergm` (see details).

- set.control.ergm.ego:

  Control arguments passed to `ergm.ego` (see details).

- verbose:

  If `TRUE`, print model fitting progress to console.

- nested.edapprox:

  Logical. If `edapprox = TRUE` the dissolution model is an initial
  segment of the formation model (see details).

- ...:

  Additional arguments passed to other functions.

## Value

A fitted network model object of class `netest`.

## Details

`netest` is a wrapper function for the `ergm`, `ergm.ego`, and `tergm`
functions that estimate static and dynamic network models. Network model
estimation is the first step in simulating a stochastic network epidemic
model in `EpiModel`. The output from `netest` is a necessary input for
running the epidemic simulations in
[`netsim`](http://epimodel.github.io/EpiModel/reference/netsim.md). With
a fitted network model, one should always first proceed to model
diagnostics, available through the
[`netdx`](http://epimodel.github.io/EpiModel/reference/netdx.md)
function, to check model fit. A detailed description of fitting these
models, along with examples, may be found in the [Network Modeling for
Epidemics](https://epimodel.github.io/sismid/) tutorials.

## Edges Dissolution Approximation

The edges dissolution approximation method is described in Carnegie et
al. This approximation requires that the dissolution coefficients are
known, that the formation model is being fit to cross-sectional data
conditional on those dissolution coefficients, and that the terms in the
dissolution model are a subset of those in the formation model. Under
certain additional conditions, the formation coefficients of a STERGM
model are approximately equal to the coefficients of that same model fit
to the observed cross-sectional data as an ERGM, minus the corresponding
coefficients in the dissolution model. The approximation thus estimates
this ERGM (which is typically much faster than estimating a STERGM) and
subtracts the dissolution coefficients.

The conditions under which this approximation best hold are when there
are few relational changes from one time step to another; i.e. when
either average relational durations are long, or density is low, or
both. Conveniently, these are the same conditions under which STERGM
estimation is slowest. Note that the same approximation is also used to
obtain starting values for the STERGM estimate when the latter is being
conducted. The estimation does not allow for calculation of standard
errors, p-values, or likelihood for the formation model; thus, this
approach is of most use when the main goal of estimation is to drive
dynamic network simulations rather than to conduct inference on the
formation model. The user is strongly encouraged to examine the behavior
of the resulting simulations to confirm that the approximation is
adequate for their purposes. For an example, see the vignette for the
package `tergm`.

It has recently been found that subtracting a modified version of the
dissolution coefficients from the formation coefficients provides a more
principled approximation, and this is now the form of the approximation
applied by `netest`. The modified values subtracted from the formation
coefficients are equivalent to the (crude) dissolution coefficients with
their target durations increased by 1. The `nested.edapprox` argument
toggles whether to implement this modified version by appending the
dissolution terms to the formation model and appending the relevant
values to the vector of formation model coefficients (value = `FALSE`),
whereas the standard version subtracts the relevant values from the
initial formation model coefficients (value = `TRUE`).

## Control Arguments

The `ergm`, `ergm.ego`, and `tergm` functions allow control settings for
the model fitting process. When fitting a STERGM directly (setting
`edapprox` to `FALSE`), control parameters may be passed to the `tergm`
function with the `set.control.tergm` argument in `netest`. The controls
should be input through the
[`control.tergm()`](https://rdrr.io/pkg/tergm/man/control.tergm.html)
function, with the available parameters listed in the
[`tergm::control.tergm`](https://rdrr.io/pkg/tergm/man/control.tergm.html)
help page in the `tergm` package.

When fitting a STERGM indirectly (setting `edapprox` to `TRUE`), control
settings may be passed to the `ergm` function using `set.control.ergm`,
or to the `ergm.ego` function using `set.control.ergm.ego`. The controls
should be input through the `control.ergm()` and `control.ergm.ego()`
functions, respectively, with the available parameters listed in the
[`ergm::control.ergm`](https://rdrr.io/pkg/ergm/man/control.ergm.html)
help page in the `ergm` package and the
[`ergm.ego::control.ergm.ego`](https://rdrr.io/pkg/ergm.ego/man/control.ergm.ego.html)
help page in the `ergm.ego` package. An example is below.

## References

Krivitsky PN, Handcock MS. "A Separable Model for Dynamic Networks."
JRSS(B). 2014; 76.1: 29-46.

Carnegie NB, Krivitsky PN, Hunter DR, Goodreau SM. An Approximation
Method for Improving Dynamic Network Model Fitting. Journal of
Computational and Graphical Statistics. 2014; 24(2): 502-519.

Jenness SM, Goodreau SM and Morris M. EpiModel: An R Package for
Mathematical Modeling of Infectious Disease over Networks. Journal of
Statistical Software. 2018; 84(8): 1-47.

## See also

Use [`netdx`](http://epimodel.github.io/EpiModel/reference/netdx.md) to
diagnose the fitted network model, and
[`netsim`](http://epimodel.github.io/EpiModel/reference/netsim.md) to
simulate epidemic spread over a simulated dynamic network consistent
with the model fit.

## Examples

``` r
# Initialize a network of 100 nodes
nw <- network_initialize(n = 100)

# Set formation formula
formation <- ~edges + concurrent

# Set target statistics for formation
target.stats <- c(50, 25)

# Obtain the offset coefficients
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)

# Estimate the STERGM using the edges dissolution approximation
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
#> Warning: ‘glpk’ selected as the solver, but package ‘Rglpk’ is not available; falling back to ‘lpSolveAPI’. This should be fine unless the sample size and/or the number of parameters is very big.
#> 1 
#> Optimizing with step length 1.0000.
#> The log-likelihood improved by 0.0835.
#> Convergence test p-value: 0.0002. 
#> Converged with 99% confidence.
#> Finished MCMLE.
#> This model was fit using MCMC.  To examine model diagnostics and check
#> for degeneracy, use the mcmc.diagnostics() function.
est
#> EpiModel Network Estimation
#> =======================
#> Model class: netest
#> Estimation Method: ERGM with Edges Approximation
#> 
#> Model Form
#> -----------------------
#> Formation: ~edges + concurrent
#> <environment: 0x55fba1f79198>
#> Target Statistics: 50 25
#> Constraints: ~.
#> 
#> Dissolution: ~offset(edges)
#> Target Statistics: 10

# To estimate the STERGM directly, use edapprox = FALSE
# est2 <- netest(nw, formation, target.stats, coef.diss, edapprox = FALSE)
```
