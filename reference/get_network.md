# Extract Network Objects from Network Simulations

Extracts the network object from either a network epidemic model object
generated with `netsim`, a network diagnostic simulation generated with
`netdx`, or a `netsim_dat` object used internally in `netsim`. For
`netdx` or `netsim` with `tergmLite == FALSE`, the extracted network
object is a `networkDynamic`, which can be collapsed down to a static
`network` object with the `collapse` and `at` arguments. For `netsim`
with `tergmLite == TRUE`, the extracted network object is the final
`networkLite`, the `collapse` argument should be `FALSE`, and the `at`
argument should be missing. For `netsim_dat`, the `collapse` and `at`
arguments are not supported, and the network object is either the
current `networkLite` (if `tergmLite == TRUE`) or the current
`networkDynamic` (if `tergmLite == FALSE`).

## Usage

``` r
get_network(x, ...)

# S3 method for class 'netdx'
get_network(x, sim = 1, collapse = FALSE, at = NULL, ...)

# S3 method for class 'netsim'
get_network(x, sim = 1, network = 1, collapse = FALSE, at = NULL, ...)

# S3 method for class 'netsim_dat'
get_network(x, network = 1L, ...)
```

## Arguments

- x:

  An `EpiModel` object of class
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md),
  [`netdx()`](http://epimodel.github.io/EpiModel/reference/netdx.md), or
  `netsim_dat`.

- ...:

  Additional arguments.

- sim:

  Simulation number of extracted network, for `netdx` and `netsim`.

- collapse:

  If `TRUE`, collapse the `networkDynamic` object to a static `network`
  object at a specified time step. Applicable to `netdx` objects and
  `netsim` objects with `tergmLite == FALSE`.

- at:

  If `collapse` is `TRUE`, the time step at which the extracted network
  should be collapsed. Applicable to `netdx` objects and `netsim`
  objects with `tergmLite == FALSE`.

- network:

  Network number, for `netsim` or `netsim_dat` objects with multiple
  overlapping networks (advanced use, and not applicable to `netdx`
  objects).

## Value

For `netdx` or `netsim` with `tergmLite == FALSE`, a `networkDynamic`
object (if `collapse = FALSE`) or a static `network` object (if
`collapse = TRUE`). For `netsim` with `tergmLite == TRUE` or
`netsim_dat` with `tergmLite == TRUE`, a `networkLite` object. For
`netsim_dat` with `tergmLite == FALSE`, a `networkDynamic` object.

## Details

This function requires that the network object is saved during the
network simulation while running either
[`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md) or
[`netdx()`](http://epimodel.github.io/EpiModel/reference/netdx.md). For
the former, that is specified by setting the `save.network` parameter in
[`control.net()`](http://epimodel.github.io/EpiModel/reference/control.net.md)
to `TRUE`. For the latter, that is specified with the `keep.tnetwork`
parameter directly in
[`netdx()`](http://epimodel.github.io/EpiModel/reference/netdx.md).

## Examples

``` r
# Set up network and TERGM formula
nw <- network_initialize(n = 100)
nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)

# Estimate the model
est <- netest(nw, formation, target.stats, coef.diss)
#> Starting maximum pseudolikelihood estimation (MPLE):
#> Obtaining the responsible dyads.
#> Evaluating the predictor and response matrix.
#> Maximizing the pseudolikelihood.
#> Finished MPLE.

# Run diagnostics, saving the networkDynamic objects
dx <- netdx(est, nsteps = 10, nsims = 3, keep.tnetwork = TRUE,
            verbose = FALSE)
#> Warning: NAs introduced by coercion to integer range
#> Warning: NAs introduced by coercion to integer range
#> Warning: NAs introduced by coercion to integer range

# Extract the network for simulation 2 from dx object
get_network(dx, sim = 2)
#> NetworkDynamic properties:
#>   distinct change times: 12 
#>   maximal time range: 0 until  Inf 
#> 
#> Includes optional net.obs.period attribute:
#>  Network observation period info:
#>   Number of observation spells: 2 
#>   Maximal time range observed: 0 until 11 
#>   Temporal mode: discrete 
#>   Time unit: step 
#>   Suggested time increment: 1 
#> 
#>  Network attributes:
#>   vertices = 100 
#>   directed = FALSE 
#>   hyper = FALSE 
#>   loops = FALSE 
#>   multiple = FALSE 
#>   bipartite = FALSE 
#>   net.obs.period: (not shown)
#>   total edges= 75 
#>     missing edges= 0 
#>     non-missing edges= 75 
#> 
#>  Vertex attribute names: 
#>     group vertex.names 
#> 
#>  Edge attribute names: 
#>     active 

# Extract and collapse the network from simulation 1 at time step 5
get_network(dx, collapse = TRUE, at = 5)
#>  Network attributes:
#>   vertices = 100 
#>   directed = FALSE 
#>   hyper = FALSE 
#>   loops = FALSE 
#>   multiple = FALSE 
#>   bipartite = FALSE 
#>   total edges= 55 
#>     missing edges= 0 
#>     non-missing edges= 55 
#> 
#>  Vertex attribute names: 
#>     group vertex.names 
#> 
#> No edge attributes

# Parameterize the epidemic model, and simulate it
param <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.15)
init <- init.net(i.num = 10, i.num.g2 = 10)
control <- control.net(type = "SI", nsteps = 10, nsims = 3, verbose = FALSE)
mod <- netsim(est, param, init, control)

# Extract the network for simulation 2 from mod object
get_network(mod, sim = 2)
#> NetworkDynamic properties:
#>   distinct change times: 12 
#>   maximal time range: 0 until  Inf 
#> 
#>  Dynamic (TEA) attributes:
#>   Vertex TEAs:    testatus.active 
#> 
#> Includes optional net.obs.period attribute:
#>  Network observation period info:
#>   Number of observation spells: 2 
#>   Maximal time range observed: 0 until 11 
#>   Temporal mode: discrete 
#>   Time unit: step 
#>   Suggested time increment: 1 
#> 
#>  Network attributes:
#>   vertices = 100 
#>   directed = FALSE 
#>   hyper = FALSE 
#>   loops = FALSE 
#>   multiple = FALSE 
#>   bipartite = FALSE 
#>   net.obs.period: (not shown)
#>   vertex.pid = tergm_pid 
#>   total edges= 81 
#>     missing edges= 0 
#>     non-missing edges= 81 
#> 
#>  Vertex attribute names: 
#>     active group status tergm_pid testatus.active vertex.names 
#> 
#>  Edge attribute names: 
#>     active 

## Extract and collapse the network from simulation 1 at time step 5
get_network(mod, collapse = TRUE, at = 5)
#>  Network attributes:
#>   vertices = 100 
#>   directed = FALSE 
#>   hyper = FALSE 
#>   loops = FALSE 
#>   multiple = FALSE 
#>   bipartite = FALSE 
#>   vertex.pid = tergm_pid 
#>   total edges= 49 
#>     missing edges= 0 
#>     non-missing edges= 49 
#> 
#>  Vertex attribute names: 
#>     group status tergm_pid testatus vertex.names 
#> 
#> No edge attributes
```
