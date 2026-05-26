# Get a Cumulative Edgelist From a Specified Network

Returns the cumulative edgelist for a given network: a record of every
edge EpiModel has tracked, with the time steps at which each edge formed
and (if no longer active) dissolved. This is the canonical mechanism for
querying partnership histories during or after a simulation, and is
particularly important under `tergmLite = TRUE`, where the full
`networkDynamic` history is not retained.

## Usage

``` r
get_cumulative_edgelist(dat, network)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md).

- network:

  Numerical index of the network from which the cumulative edgelist
  should be extracted. (May be \> 1 for models with multiple overlapping
  networks.)

## Value

A `tibble` with four columns:

- `head`: the unique ID (see
  [`get_unique_ids()`](https://epimodel.github.io/EpiModel/reference/unique_id-tools.md))
  of the head node.

- `tail`: the unique ID of the tail node.

- `start`: the time step at which the edge formed.

- `stop`: the time step at which the edge dissolved, or `NA` if the edge
  is still active. Edges are active over `[start, stop]` inclusive.

## Details

Cumulative-edgelist tracking is opt-in. It must be enabled via
`control.net(cumulative.edgelist = TRUE)` (see
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md));
calling this function on a simulation that did not enable tracking
raises an error.

Inside a custom module (where `dat` is the live `netsim_dat` object),
call this function directly. After a
[`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md)
run, the saved edgelists are attached to the returned object as
`sim$cumulative.edgelist[[s]]` (one element per simulation, populated
when `save.cumulative.edgelist = TRUE`); read those directly rather than
calling this function on the processed output.

For the *current* (active-only, no history) edgelist, see
[`get_edgelist()`](https://epimodel.github.io/EpiModel/reference/get_edgelist.md)
and
[`get_edgelists_df()`](https://epimodel.github.io/EpiModel/reference/get_edgelists_df.md).

## See also

[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md)
for the controlling flags (`cumulative.edgelist`, `truncate.el.cuml`,
`save.cumulative.edgelist`).
[`vignette("network-objects", package = "EpiModel")`](https://epimodel.github.io/EpiModel/articles/network-objects.md)
walks through the full lifecycle.

Other cumulative_edgelist:
[`as_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/as_cumulative_edgelist.md),
[`dedup_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/dedup_cumulative_edgelist.md),
[`get_cumulative_degree()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_degree.md),
[`get_cumulative_edgelists_df()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_edgelists_df.md),
[`get_partners()`](https://epimodel.github.io/EpiModel/reference/get_partners.md),
[`reachable-nodes`](https://epimodel.github.io/EpiModel/reference/reachable-nodes.md),
[`update_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/update_cumulative_edgelist.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Inside a custom module (dat is the live netsim_dat object):
el <- get_cumulative_edgelist(dat, network = 1)

# Post-simulation on a processed netsim object, read the saved slot:
sim$cumulative.edgelist[[1]]
} # }
```
