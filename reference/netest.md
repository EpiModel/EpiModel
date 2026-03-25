# Dynamic Network Model Estimation

Estimates statistical network models using the exponential random graph
modeling (ERGM) framework with extensions for dynamic/temporal models
(STERGM). This is typically the first step in the network modeling
pipeline, followed by
[`netdx`](https://epimodel.github.io/EpiModel/reference/netdx.md) for
model diagnostics and
[`netsim`](https://epimodel.github.io/EpiModel/reference/netsim.md) for
epidemic simulation.

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
  set.control.tergm = control.tergm(MCMC.maxchanges = .Machine$integer.max),
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
  where `...` are additional network statistics (ERGM terms). This
  formula specifies which structural features of the network should be
  reproduced by the model. Common terms include `edges` (overall
  connectivity), `nodematch` (homophily by attribute), `concurrent`
  (overlapping partnerships), and `degree` (degree distribution
  constraints). See
  [`ergm::ergm-terms`](https://rdrr.io/pkg/ergm/man/ergmTerm.html) for
  the full list of available terms.

- target.stats:

  Vector of target statistics for the formation model, with one number
  for each network statistic in the model. These are the observed (or
  desired) values for each term in the formation formula. For example,
  if `formation = ~edges + concurrent`, then `target.stats = c(175, 40)`
  means the model should produce approximately 175 edges and 40 nodes
  with 2 or more partners. For an `edges`-only model, a useful starting
  value is `mean_degree * network_size / 2`. Ignored if fitting via
  `ergm.ego`.

- coef.diss:

  An object of class `disscoef` output from the
  [`dissolution_coefs`](https://epimodel.github.io/EpiModel/reference/dissolution_coefs.md)
  function. This encodes the average partnership duration(s) and the
  corresponding dissolution model coefficients. For models with vital
  dynamics (arrivals and departures), the `d.rate` argument in
  [`dissolution_coefs`](https://epimodel.github.io/EpiModel/reference/dissolution_coefs.md)
  should be set to adjust for the competing risk of node departure.

- constraints:

  Right-hand sided formula specifying constraints for the modeled
  network, in the form `~...`, where `...` are constraint terms. By
  default, no constraints are set.

- coef.form:

  Vector of coefficients for the offset terms in the formation formula.

- edapprox:

  If `TRUE`, use the indirect edges dissolution approximation method for
  the dynamic model fit, otherwise use the more time-intensive full
  STERGM estimation (see details). The approximation is recommended for
  most use cases, especially when average partnership durations are
  moderate to long (\> 25 time steps). Direct STERGM estimation
  (`edapprox = FALSE`) is slower but may be preferred for very short
  durations or when inferential quantities (standard errors, p-values)
  on formation coefficients are needed. For `nw` of class `egor`, only
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

A fitted network model object of class `netest`. This object is passed
to [`netdx`](https://epimodel.github.io/EpiModel/reference/netdx.md) for
diagnostics and to
[`netsim`](https://epimodel.github.io/EpiModel/reference/netsim.md) for
epidemic simulation. Use [`print()`](https://rdrr.io/r/base/print.html)
to view the model form, including the formation formula, target
statistics, and dissolution model. Use
[`summary()`](https://rdrr.io/r/base/summary.html) to view the estimated
model coefficients and goodness-of-fit statistics from the underlying
ERGM or STERGM fit.

## Details

`netest` is a wrapper function for the `ergm`, `ergm.ego`, and `tergm`
functions that estimate static and dynamic network models. Network model
estimation is the first step in simulating a stochastic network epidemic
model in `EpiModel`. The output from `netest` is a necessary input for
running the epidemic simulations in
[`netsim`](https://epimodel.github.io/EpiModel/reference/netsim.md).
With a fitted network model, one should always first proceed to model
diagnostics, available through the
[`netdx`](https://epimodel.github.io/EpiModel/reference/netdx.md)
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

## Typical Workflow

The network modeling pipeline in EpiModel typically follows these steps:

1.  Initialize a network:
    [`network_initialize`](https://epimodel.github.io/EpiModel/reference/network_initialize.md)

2.  Specify formation and dissolution: a formation formula (e.g.,
    `~edges + concurrent`) with target statistics, and dissolution
    coefficients via
    [`dissolution_coefs`](https://epimodel.github.io/EpiModel/reference/dissolution_coefs.md)

3.  Estimate the network model: `netest()`

4.  Diagnose model fit:
    [`netdx`](https://epimodel.github.io/EpiModel/reference/netdx.md)

5.  Simulate the epidemic:
    [`netsim`](https://epimodel.github.io/EpiModel/reference/netsim.md)
    with
    [`param.net`](https://epimodel.github.io/EpiModel/reference/param.net.md),
    [`init.net`](https://epimodel.github.io/EpiModel/reference/init.net.md),
    and
    [`control.net`](https://epimodel.github.io/EpiModel/reference/control.net.md)

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

Use
[`dissolution_coefs`](https://epimodel.github.io/EpiModel/reference/dissolution_coefs.md)
to compute dissolution coefficients before estimation. Use
[`netdx`](https://epimodel.github.io/EpiModel/reference/netdx.md) to
diagnose the fitted network model, and
[`netsim`](https://epimodel.github.io/EpiModel/reference/netsim.md) to
simulate epidemic spread over a simulated dynamic network consistent
with the model fit. Parameterize the epidemic simulation with
[`param.net`](https://epimodel.github.io/EpiModel/reference/param.net.md),
[`init.net`](https://epimodel.github.io/EpiModel/reference/init.net.md),
and
[`control.net`](https://epimodel.github.io/EpiModel/reference/control.net.md).

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
#> Starting simulated annealing (SAN)
#> Iteration 1 of at most 4
#> Finished simulated annealing
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
#> The log-likelihood improved by 0.1617.
#> Estimating equations are not within tolerance region.
#> Iteration 2 of at most 60:
#> 1 
#> Optimizing with step length 1.0000.
#> The log-likelihood improved by 0.0348.
#> Convergence test p-value: 0.0001. 
#> Converged with 99% confidence.
#> Finished MCMLE.
#> This model was fit using MCMC.  To examine model diagnostics and check
#> for degeneracy, use the mcmc.diagnostics() function.
# View the model form (formation, targets, dissolution)
est
#> EpiModel Network Estimation
#> =======================
#> Model class: netest
#> Estimation Method: ERGM with Edges Approximation
#> 
#> Model Form
#> -----------------------
#> Formation: ~edges + concurrent
#> <environment: 0x55ad3f6e6e98>
#> Target Statistics: 50 25
#> Constraints: ~.
#> 
#> Dissolution: ~offset(edges)
#> Target Statistics: 10

# View the estimated coefficients
summary(est)
#> Call:
#> ergm(formula = formation, constraints = constraints, offset.coef = coef.form, 
#>     target.stats = target.stats, eval.loglik = FALSE, control = set.control.ergm, 
#>     verbose = verbose, basis = nw)
#> 
#> Monte Carlo Maximum Likelihood Results:
#> 
#>            Estimate Std. Error MCMC % z value Pr(>|z|)    
#> edges       -4.4092     0.3236      0 -13.627   <1e-04 ***
#> concurrent  -0.2635     0.4044      0  -0.652    0.515    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood was not estimated for this fit. To get deviances, AIC, and/or BIC, use ‘*fit* <-logLik(*fit*, add=TRUE)’ to add it to the object or rerun this function with eval.loglik=TRUE.
#> 
#> Dissolution Coefficients
#> =======================
#> Dissolution Model: ~offset(edges)
#> Target Statistics: 10
#> Crude Coefficient: 2.197225
#> Mortality/Exit Rate: 0
#> Adjusted Coefficient: 2.197225

if (FALSE) { # \dontrun{
# Model with homophily on a nodal attribute
nw2 <- network_initialize(n = 500)
nw2 <- set_vertex_attribute(nw2, "risk", rep(0:1, each = 250))
formation2 <- ~edges + nodematch("risk")
target.stats2 <- c(175, 110)
coef.diss2 <- dissolution_coefs(dissolution = ~offset(edges), duration = 50)
est2 <- netest(nw2, formation2, target.stats2, coef.diss2)
est2

# Direct STERGM estimation (slower, for short durations or inference)
est3 <- netest(nw, formation, target.stats, coef.diss, edapprox = FALSE)
} # }
```
