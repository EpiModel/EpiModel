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
  It may be omitted when the model has exactly one clique layer.

## Value

The updated `netsim_dat` object.

## Details

The clique edges are built from the grouping attribute once, when the
simulation starts, and are afterwards changed only by arrivals,
departures, and this function. Setting the grouping attribute with
[`set_attr()`](https://epimodel.github.io/EpiModel/reference/net-accessor.md)
alone therefore does not rewire the layer: the old edges remain and the
layer no longer matches the groups. `move_to_group()` sets the attribute
and the edges together, in either network storage mode
(`tergmLite = TRUE` or `FALSE`), so that the layer stays exactly the
union of the cliques of the current groups.

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
if (FALSE) { # \dontrun{
# A module in which people aged 18 to 29 leave home with a monthly
# probability of 0.01, each starting a household of their own
leave_home <- function(dat, at) {
  active <- get_attr(dat, "active")
  age <- get_attr(dat, "age")
  hh_id <- get_attr(dat, "hh_id")
  elig <- which(active == 1 & age >= 18 & age < 30)
  movers <- elig[runif(length(elig)) < 0.01]
  if (length(movers) > 0) {
    new_ids <- max(hh_id, na.rm = TRUE) + seq_along(movers)
    dat <- move_to_group(dat, ids = movers, group = new_ids)
  }
  return(dat)
}
} # }
```
