# Get the Cumulative Edgelists of a Model

Combines the cumulative edgelists from one or more network layers into a
single `data.frame`, adding a `network` column to identify the layer.
Like
[`get_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_edgelist.md),
this requires `control.net(cumulative.edgelist = TRUE)`.

## Usage

``` r
get_cumulative_edgelists_df(dat, networks = NULL)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md).

- networks:

  Numerical indexes of the networks to extract the partnerships from.
  (May be \> 1 for models with multiple overlapping networks.) If
  `NULL`, extract from all networks.

## Value

A `data.frame` with 5 columns:

- `head`: the unique ID (see
  [`get_unique_ids()`](https://epimodel.github.io/EpiModel/reference/unique_id-tools.md))
  of the head node.

- `tail`: the unique ID of the tail node.

- `start`: the time step at which the edge formed.

- `stop`: the time step at which the edge dissolved, or `NA` if still
  active.

- `network`: the numerical index of the network on which the partnership
  lives.

Note: column names `head`/`tail` here match the single-network
[`get_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_edgelist.md)
output. The `index`/`partner` naming used by
[`get_partners()`](https://epimodel.github.io/EpiModel/reference/get_partners.md)
is reserved for queries about specific index nodes.

## See also

[`get_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_edgelist.md)
(single network),
[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md)
(`cumulative.edgelist`, `save.cumulative.edgelist`).
[`vignette("network-objects", package = "EpiModel")`](https://epimodel.github.io/EpiModel/articles/network-objects.md)
walks through the full lifecycle.

Other cumulative_edgelist:
[`as_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/as_cumulative_edgelist.md),
[`dedup_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/dedup_cumulative_edgelist.md),
[`get_cumulative_degree()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_degree.md),
[`get_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_edgelist.md),
[`get_partners()`](https://epimodel.github.io/EpiModel/reference/get_partners.md),
[`reachable-nodes`](https://epimodel.github.io/EpiModel/reference/reachable-nodes.md),
[`update_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/update_cumulative_edgelist.md)
