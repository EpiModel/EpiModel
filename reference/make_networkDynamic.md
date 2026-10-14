# Build a networkDynamic Object from Simulation Output

Reconstructs a `networkDynamic` object from a completed `netsim`
simulation, using the cumulative edgelist for edge spells and the
tracked attribute history for vertex activity and time-varying
attributes. The resulting object can be used with `ndtv` for
visualization or with `tsna` for temporal network analysis.

## Usage

``` r
make_networkDynamic(sim, sim_num = 1, network = NULL)
```

## Arguments

- sim:

  An `EpiModel` object of class `netsim`.

- sim_num:

  The simulation number to use (default 1).

- network:

  Optional network number to filter edges. If `NULL` (default), edges
  from all networks are included.

## Value

A `networkDynamic` object with active vertex spells and time-varying
vertex attributes.

## Details

The simulation must be run with `cumulative.edgelist = TRUE`,
`save.cumulative.edgelist = TRUE`, and `truncate.el.cuml = Inf` in
[`control.net`](https://epimodel.github.io/EpiModel/reference/control.net.md),
so that the cumulative edgelist holds every edge, and with at least one
attribute in `tracked.attributes` or `tracked.attributes.once`, so that
the activity of the nodes is recorded. The vertex ids of the
`networkDynamic` object are the `unique_id`s of the nodes. Each pair of
nodes tied on a layer has one edge, with one activity spell per period
the tie was active, and the layer in the `network` edge attribute. With
`network = NULL` and more than one layer, the layers are merged into one
network that allows multiple edges (`multiple = TRUE`), so that a pair
of nodes tied on two layers has one edge per layer.

The nodes are active from the first record of the attribute history: for
a new simulation from time 0, as in the `networkDynamic` objects that
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
saves without `tergmLite`, and for a restarted one from `start - 1`, the
snapshot its history starts from. An edge that formed before the
activity spell of one of its nodes starts is active from the start of
that spell, so that no edge is active between inactive vertices.

## See also

[`get_attr_at()`](https://epimodel.github.io/EpiModel/reference/get_attr_at.md),
[`get_attr_history()`](https://epimodel.github.io/EpiModel/reference/get_attr_history.md)

## Examples

``` r
# \donttest{
nw <- network_initialize(n = 100)
est <- netest(nw, formation = ~edges, target.stats = 50,
              coef.diss = dissolution_coefs(~offset(edges), 10),
              verbose = FALSE)
#> Starting simulated annealing (SAN)
#> Iteration 1 of at most 4
#> Finished simulated annealing
#> Starting maximum pseudolikelihood estimation (MPLE):
#> Obtaining the responsible dyads.
#> Evaluating the predictor and response matrix.
#> Maximizing the pseudolikelihood.
#> Finished MPLE.
control <- control.net(type = "SI", nsteps = 20, nsims = 1,
                       tergmLite = TRUE, resimulate.network = TRUE,
                       cumulative.edgelist = TRUE,
                       save.cumulative.edgelist = TRUE,
                       truncate.el.cuml = Inf,
                       tracked.attributes = "status", verbose = FALSE)
sim <- netsim(est, param.net(inf.prob = 0.3, act.rate = 1),
              init.net(i.num = 10), control)

# The full network history of a tergmLite run, with each node's status
nd <- make_networkDynamic(sim)
nd
#> NetworkDynamic properties:
#>   distinct change times: 22 
#>   maximal time range: 0 until  Inf 
#> 
#>  Dynamic (TEA) attributes:
#>   Vertex TEAs:    status.active 
#> 
#>  Network attributes:
#>   vertices = 100 
#>   directed = FALSE 
#>   hyper = FALSE 
#>   loops = FALSE 
#>   multiple = FALSE 
#>   bipartite = FALSE 
#>   total edges= 165 
#>     missing edges= 0 
#>     non-missing edges= 165 
#> 
#>  Vertex attribute names: 
#>     active status.active vertex.names 
#> 
#>  Edge attribute names: 
#>     active network 
x10 <- network.collapse(nd, at = 10)
table(get_vertex_attribute(x10, "status"))
#> 
#>  i  s 
#> 40 60 
# }
```
