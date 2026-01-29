# Dynamic Network Model Diagnostics

Runs dynamic diagnostics on an ERGM/STERGM estimated with
[`netest`](http://epimodel.github.io/EpiModel/reference/netest.md).

## Usage

``` r
netdx(
  x,
  nsims = 1,
  dynamic = TRUE,
  nsteps = NULL,
  nwstats.formula = "formation",
  set.control.ergm = control.simulate.formula(),
  set.control.tergm = control.simulate.formula.tergm(MCMC.maxchanges = Inf),
  sequential = TRUE,
  keep.tedgelist = FALSE,
  keep.tnetwork = FALSE,
  verbose = TRUE,
  ncores = 1,
  skip.dissolution = FALSE
)
```

## Arguments

- x:

  An `EpiModel` object of class `netest`.

- nsims:

  Number of simulations to run.

- dynamic:

  If `TRUE`, runs dynamic diagnostics. If `FALSE` and the `netest`
  object was fit with the Edges Dissolution approximation method,
  simulates from the static ERGM fit.

- nsteps:

  Number of time steps per simulation (dynamic simulations only).

- nwstats.formula:

  A right-hand sided ERGM formula with the network statistics of
  interest. The default is the formation formula of the network model
  contained in `x`.

- set.control.ergm:

  Control arguments passed to `ergm`'s `simulate_formula.network` (see
  details).

- set.control.tergm:

  Control arguments passed to `tergm`'s `simulate_formula.network` (see
  details).

- sequential:

  For static diagnostics (`dynamic=FALSE`): if `FALSE`, each of the
  `nsims` simulated Markov chains begins at the initial network; if
  `TRUE`, the end of one simulation is used as the start of the next.

- keep.tedgelist:

  If `TRUE`, keep the timed edgelist generated from the dynamic
  simulations. Returned in the form of a list of matrices, with one
  entry per simulation. Accessible at `$edgelist`.

- keep.tnetwork:

  If `TRUE`, keep the full networkDynamic objects from the dynamic
  simulations. Returned in the form of a list of nD objects, with one
  entry per simulation. Accessible at `$network`.

- verbose:

  If `TRUE`, print progress to the console.

- ncores:

  Number of processor cores to run multiple simulations on, using the
  `foreach` and `doParallel` implementations.

- skip.dissolution:

  If `TRUE`, skip over the calculations of duration and dissolution
  stats in `netdx`.

## Value

A list of class `netdx`.

## Details

The `netdx` function handles dynamic network diagnostics for network
models fit with the
[`netest`](http://epimodel.github.io/EpiModel/reference/netest.md)
function. Given the fitted model, `netdx` simulates a specified number
of dynamic networks for a specified number of time steps per simulation.
The network statistics in `nwstats.formula` are saved for each time
step. Summary statistics for the formation model terms, as well as
dissolution model and relational duration statistics, are then
calculated and can be accessed when printing or plotting the `netdx`
object. See
[`print.netdx`](http://epimodel.github.io/EpiModel/reference/print.netdx.md)
and
[`plot.netdx`](http://epimodel.github.io/EpiModel/reference/plot.netdx.md)
for details on printing and plotting.

## Control Arguments

Models fit with the full STERGM method in `netest` (setting the
`edapprox` argument to `FALSE`) require only a call to `tergm`'s
`simulate_formula.network`. Control parameters for those simulations may
be set using `set.control.tergm` in `netdx`. The parameters should be
input through the `control.simulate.formula.tergm` function, with the
available parameters listed in the
[`tergm::control.simulate.formula.tergm`](https://rdrr.io/pkg/tergm/man/control.simulate.tergm.html)
help page in the `tergm` package.

Models fit with the ERGM method with the edges dissolution approximation
(setting `edapprox` to `TRUE`) require a call first to `ergm`'s
`simulate_formula.network` for simulating an initial network, and second
to `tergm`'s `simulate_formula.network` for simulating that static
network forward through time. Control parameters may be set for both
processes in `netdx`. For the first, the parameters should be input
through the `control.simulate.formula()` function, with the available
parameters listed in the
[`ergm::control.simulate.formula`](https://rdrr.io/pkg/ergm/man/control.simulate.ergm.html)
help page in the `ergm` package. For the second, parameters should be
input through the
[`control.simulate.formula.tergm()`](https://rdrr.io/pkg/tergm/man/control.simulate.tergm.html)
function, with the available parameters listed in the
[`tergm::control.simulate.formula.tergm`](https://rdrr.io/pkg/tergm/man/control.simulate.tergm.html)
help page in the `tergm` package. An example is shown below.

## See also

Plot these model diagnostics with
[`plot.netdx`](http://epimodel.github.io/EpiModel/reference/plot.netdx.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Network initialization and model parameterization
nw <- network_initialize(n = 100)
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 25)

# Estimate the model
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

# Static diagnostics on the ERGM fit
dx1 <- netdx(est,
  nsims = 1e4, dynamic = FALSE,
  nwstats.formula = ~ edges + meandeg + concurrent
)
dx1
plot(dx1, method = "b", stats = c("edges", "concurrent"))

# Dynamic diagnostics on the STERGM approximation
dx2 <- netdx(est,
  nsims = 5, nsteps = 500,
  nwstats.formula = ~ edges + meandeg + concurrent,
  set.control.ergm = control.simulate.formula(MCMC.burnin = 1e6)
)
dx2
plot(dx2, stats = c("edges", "meandeg"), plots.joined = FALSE)
plot(dx2, type = "duration")
plot(dx2, type = "dissolution", qnts.col = "orange2")
plot(dx2, type = "dissolution", method = "b", col = "bisque")

# Dynamic diagnostics on a more complex model
nw <- network_initialize(n = 1000)
nw <- set_vertex_attribute(nw, "neighborhood", rep(1:10, 100))
formation <- ~edges + nodematch("neighborhood", diff = TRUE)
target.stats <- c(800, 45, 81, 24, 16, 32, 19, 42, 21, 24, 31)
coef.diss <- dissolution_coefs(dissolution = ~offset(edges) +
                    offset(nodematch("neighborhood", diff = TRUE)),
                    duration = c(52, 58, 61, 55, 81, 62, 52, 64, 52, 68, 58))
est2 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
dx3 <- netdx(est2, nsims = 5, nsteps = 100)
print(dx3)
plot(dx3)
plot(dx3, type = "duration", plots.joined = TRUE, qnts = 0.2, legend = TRUE)
plot(dx3, type = "dissolution", mean.smooth = FALSE, mean.col = "red")
} # }
```
