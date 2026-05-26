# Return the Cumulative Degree of a Set of Index Nodes

Counts the number of distinct partners each of a set of index nodes has
accumulated over the tracked history. A thin wrapper around
[`get_partners()`](https://epimodel.github.io/EpiModel/reference/get_partners.md)
that collapses the partner list to a count per index.

## Usage

``` r
get_cumulative_degree(
  dat,
  index_posit_ids,
  networks = NULL,
  truncate = Inf,
  only.active.nodes = FALSE
)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md).

- index_posit_ids:

  The positional IDs of the indexes of interest.

- networks:

  Numerical indexes of the networks to extract the partnerships from.
  (May be \> 1 for models with multi-layer networks.) If `NULL`, extract
  from all networks.

- truncate:

  After how many time steps a partnership that is no longer active
  should be removed from the output. See the Truncation section.

- only.active.nodes:

  If `TRUE`, then inactive (e.g., deceased) partners will be removed
  from the output.

## Value

A `data.frame` with 2 columns:

- `index_pid`: the positional ID (see
  [`get_posit_ids()`](https://epimodel.github.io/EpiModel/reference/unique_id-tools.md))
  of the indexes.

- `degree`: the cumulative degree of the index.

## Cumulative Degree

The cumulative degree of a node is the number of distinct partners
connected to it during the tracked time window. The window is whatever
history the cumulative edgelist currently contains (controlled at
simulation time by `control$truncate.el.cuml` and at call time by the
`truncate` argument).

## See also

[`get_partners()`](https://epimodel.github.io/EpiModel/reference/get_partners.md)
for the full partner list backing this count.

Other cumulative_edgelist:
[`as_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/as_cumulative_edgelist.md),
[`dedup_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/dedup_cumulative_edgelist.md),
[`get_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_edgelist.md),
[`get_cumulative_edgelists_df()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_edgelists_df.md),
[`get_partners()`](https://epimodel.github.io/EpiModel/reference/get_partners.md),
[`reachable-nodes`](https://epimodel.github.io/EpiModel/reference/reachable-nodes.md),
[`update_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/update_cumulative_edgelist.md)

## Examples

``` r
if (FALSE) { # \dontrun{
get_cumulative_degree(dat, index_posit_ids = 1:50, truncate = 52)
} # }
```
