# Return the Cumulative Degree of a Set of Index Nodes

Return the Cumulative Degree of a Set of Index Nodes

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
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

- index_posit_ids:

  The positional IDs of the indexes of interest.

- networks:

  Numerical indexes of the networks to extract the partnerships from.
  (May be \> 1 for models with multi-layer networks.) If `NULL`, extract
  from all networks.

- truncate:

  After how many time steps a partnership that is no longer active
  should be removed from the output.

- only.active.nodes:

  If `TRUE`, then inactive (e.g., deceased) partners will be removed
  from the output.

## Value

A `data.frame` with 2 columns:

- `index_pid`: the positional ID (see `get_posit_ids`) of the indexes.

- `degree`: the cumulative degree of the index.

## Cumulative Degree

The cumulative degree of a node is the number of edges connected to this
node at during the time window. The time window is by default all the
steps stored in the `cumulative_edgelist` or set by the `truncate`
parameter.
