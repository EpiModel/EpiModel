# Move Nodes to Another Group of a Clique Layer

Changes the group of existing nodes during a simulation and rewires a
[`netclique()`](https://epimodel.github.io/EpiModel/reference/netclique.md)
layer to match: each moved node loses its edges to the members of its
old group and is connected to every active member of its new group.
Custom modules use it to represent people changing groups, such as a
young adult leaving home, a couple forming a household, or a child
moving in with relatives.

## Usage

``` r
move_to_group(dat, ids, group, network = NULL)
```

## Arguments

- dat:

  Main `netsim_dat` object passed through
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
  calls.

- ids:

  Positional ids of the nodes to move. Every node must be active.

- group:

  New group ids: one per element of `ids`, or a single value for all of
  them. An id in use joins that group; an id not in use starts a new
  group, which nodes given the same new id in one call form together;
  `NA` takes the node out of its group and leaves it as an isolate on
  the layer.

- network:

  Index of the clique layer in the list of layers passed to
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md).
  It may be omitted when the model has exactly one clique layer. Every
  other clique layer built on the same grouping attribute is rewired as
  well.

## Value

The updated `netsim_dat` object.

## Details

The clique edges are built from the grouping attribute once, when the
simulation starts, and are afterwards changed only by arrivals,
departures, and this function. Setting the grouping attribute with
[`set_attr()`](https://epimodel.github.io/EpiModel/reference/get_attr.md)
alone therefore does not rewire the layer: the old edges remain and the
layer no longer matches the groups. `move_to_group()` sets the attribute
and the edges together, in either network storage mode
(`tergmLite = TRUE` or `FALSE`), so that the layer stays exactly the
union of the cliques of the current groups.

Nodes that arrive in a time step join the layer when
[`nwupdate.net()`](https://epimodel.github.io/EpiModel/reference/nwupdate.net.md)
runs, which places them by the layer's arrival rule (see
[`netclique()`](https://epimodel.github.io/EpiModel/reference/netclique.md)),
so they cannot be moved before then. To place arriving nodes in chosen
groups, set the grouping attribute in the arrivals module and use
`arrivals = "join"` with an `arrivals.FUN` that returns that attribute
for the new nodes.

A node moved to its own current group keeps its edges. When the moves
leave a group without members, the group ceases to exist, as when its
last member departs. For numeric group ids, one more than the largest id
in use, `max(get_attr(dat, group.attr), na.rm = TRUE) + 1`, is an unused
id; ids of groups whose members have all departed are also unused, so
reusing one starts a new group. The ids set here count toward the
running maximum from which the `"new"` arrival rule of
[`netclique()`](https://epimodel.github.io/EpiModel/reference/netclique.md)
draws, so a later arrival is never given the id of a group formed by a
move.

Under `tergmLite = TRUE` with `tergmLite.track.duration = TRUE`, the
moved nodes' new edges are recorded as formed at the current time step.
In the cumulative edgelist (see
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md)),
a move ends the node's old edges and starts new ones.

## See also

[`netclique()`](https://epimodel.github.io/EpiModel/reference/netclique.md)
for the clique layer and its rules for placing arriving nodes.

## Examples

``` r
# \donttest{
# A module in which people aged 18 to 29 leave home with a probability of
# 0.05 per time step, each starting a household of their own
leave_home <- function(dat, at) {
  active <- get_attr(dat, "active")
  age <- get_attr(dat, "age")
  hh_id <- get_attr(dat, "hh_id")
  elig <- which(active == 1 & age >= 18 & age < 30)
  movers <- elig[runif(length(elig)) < 0.05]
  if (length(movers) > 0) {
    new_ids <- max(hh_id, na.rm = TRUE) + seq_along(movers)
    dat <- move_to_group(dat, ids = movers, group = new_ids)
  }
  return(dat)
}

# 50 households of 4 under a community layer
set.seed(12)
nw <- network_initialize(n = 200)
nw <- set_vertex_attribute(nw, "age", sample(0:79, 200, replace = TRUE))
nw <- set_vertex_attribute(nw, "hh_id", rep(1:50, each = 4))
hh <- netclique(nw, group.attr = "hh_id")
est <- netest(nw, formation = ~edges, target.stats = 100,
              coef.diss = dissolution_coefs(~offset(edges), 20),
              verbose = FALSE)
#> Starting simulated annealing (SAN)
#> Iteration 1 of at most 4
#> Finished simulated annealing
#> Starting maximum pseudolikelihood estimation (MPLE):
#> Obtaining the responsible dyads.
#> Evaluating the predictor and response matrix.
#> Maximizing the pseudolikelihood.
#> Finished MPLE.

param <- param.net(inf.prob = multilayer(0.2, 0.05), act.rate = 1)
init <- init.net(i.num = 10)
control <- control.net(type = NULL, nsteps = 20, nsims = 1,
                       infection.FUN = infection.net,
                       leave_home.FUN = leave_home,
                       tergmLite = TRUE, resimulate.network = TRUE,
                       save.run = TRUE, verbose = FALSE)
sim <- netsim(list(hh, est), param, init, control)

# More households at the end than at the start
length(unique(sim$run[[1]]$attr$hh_id))
#> [1] 65
# }
```
