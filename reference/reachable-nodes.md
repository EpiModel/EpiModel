# Get the Forward or Backward Reachable Nodes for a Set of Nodes

These functions return the Forward or Backward Reachable Nodes of a set
of nodes in a network over a time. Warning, these functions ignore nodes
without edges in the period of interest. See the `Number of Nodes`
section for details It is much faster than iterating
[`tsna::tPath`](https://rdrr.io/pkg/tsna/man/paths.html). The distance
between to each node can be back calculated using the length of the
reachable set at each time step and the fact that the reachable sets are
ordered by the time to arrival.

## Usage

``` r
get_forward_reachable(
  el_cuml,
  from_step,
  to_step,
  nodes = NULL,
  dense_optim = "auto"
)

get_backward_reachable(
  el_cuml,
  from_step,
  to_step,
  nodes = NULL,
  dense_optim = "auto"
)
```

## Arguments

- el_cuml:

  a cumulative edgelist object. That is a data.frame with at least
  columns: head, tail, start and stop. Start and stop are inclusive.

- from_step:

  the beginning of the time period.

- to_step:

  the end of the time period.

- nodes:

  the subset of nodes to calculate the FRP for. (default = NULL, all
  nodes)

- dense_optim:

  pre-process the adjacency list to speed up the computations on dense
  networks. "auto" (default), enable the optimisation when
  `n_edges > n_nodes`. "yes" always enables and "no" always disables.
  The overhead of the optimization is not worth it on sparse networks.

## Value

A named list containing: `reached`: the set of reachable nodes for each
of the `nodes`. `lengths`: A matrix of `length(nodes)` rows and one
column per timestep + 1 with the length of the reachable set at each
step from `from_step - 1` to `to_step`. The first column is always one
as the set of reachables at the beginning is just the node itself.

## Number of Nodes

To speed up the calculations and lower the memory usage, these functions
only take into account nodes with edges in the cumulative edgelist over
the period of interest. The nodes are identified in the `reached` and
`lengths` sublists by names (e.g. `node_1093`). Nodes without any edges
are therefore not calculated as the only node they reach is themselve
(length of 1). Take this fact into account when exploring the
distribution of Forward Reachable Paths for example. As the nodes with
FRP == 1 are not in the output.

## Time and Memory Use

These functions may be used to efficiently calculate multiple sets of
reachable nodes. As cumulative edgelists are way smaller than full
`networkDynamic` objects, theses functions are suited for large and
dense networks. Also, as long as the size of the `nodes` set is greater
than 5, theses functions are faster than iterating over
[`tsna::tPath`](https://rdrr.io/pkg/tsna/man/paths.html).

## Displaying Progress

These functions are using the [progressr
package](https://progressr.futureverse.org/articles/progressr-01-intro.html)
to display its progression. Use
`progressr::with_progress({ fwd_reach <- get_forward_reachable(el, from = 1, to = 260) })`
to display the progress bar. Or see the [progressr
package](https://progressr.futureverse.org/articles/progressr-01-intro.html)
for more information and customization.

## Examples

``` r
if (FALSE) { # \dontrun{

# load a network dynamic object
nd <- readRDS("nd_obj.Rds")
# convert it to a cumulative edgelist
el_cuml <- as_cumulative_edgelist(nd)

# sample 100 node indexes
nnodes <- max(el_cuml$head, el_cuml$tail)
nodes <- sample(nnodes, 100)

# `get_forward_reachable` uses steps [from_step, to_step] inclusive
el_fwd <- get_forward_reachable(el_cuml, 1, 52, nodes)[["reached"]]

# check if the results are consistent with `tsna::tPath`
nodes <- strsplit(names(el_fwd), "_")
for (i in seq_along(el_fwd)) {
  node <- as.integer(nodes[[i]][2])
  t_fwd <- tsna::tPath(
    nd, v = node,
    start = 1, end = 52 + 1, # tPath works from [start, end) right exclusive
    direction = "fwd"
  )

  t_fwd_set <- which(t_fwd$tdist < Inf)
  if(!setequal(el_fwd[[i]], t_fwd_set))
    stop("Missmatch on node: ", node)
}

# Backward:
el_bkwd <- get_backward_reachable(el_cuml, 1, 52, nodes = 1)[["reached"]]
nodes <- strsplit(names(el_bkwd), "_")
t_bkwd <- tsna::tPath(
  nd, v = nodes[i][2],
  start = 1, end = 52 + 1,
  direction = "bkwd", type = "latest.depart"
)
t_bkwd_set <- which(t_bkwd$tdist < Inf)
setequal(el_bkwd[[1]], t_bkwd_set)

} # }
```
