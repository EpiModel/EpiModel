# Return the Historical Contacts (Partners) of a Set of Index Nodes

Pulls every cumulative-edgelist row touching one of the supplied "index"
nodes, returning each `(index, partner)` pair together with the
partnership start/stop times and the network layer. This is the building
block for contact tracing over the simulated network history.

## Usage

``` r
get_partners(
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

A `data.frame` with 5 columns:

- `index`: the unique IDs of the indexes.

- `partner`: the unique IDs of the partners/contacts.

- `start`: the time step at which the edge started.

- `stop`: the time step in which the edge stopped; if ongoing, then `NA`
  is returned.

- `network`: the numerical index for the network on which the
  partnership/contact is located.

## Details

Indexes are passed as positional IDs but the output uses unique IDs,
because partners may include nodes that have already departed (and thus
no longer have a positional ID). Use
[`get_unique_ids()`](https://epimodel.github.io/EpiModel/reference/unique_id-tools.md)
and
[`get_posit_ids()`](https://epimodel.github.io/EpiModel/reference/unique_id-tools.md)
to convert between the two systems.

The `truncate` argument here filters by edge age: only edges whose
`stop` step is within `truncate` time steps of the current step are kept
(active edges are always included). It operates on whatever history is
already stored in the cumulative edgelist; it cannot recover edges that
were dropped earlier by
[`update_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/update_cumulative_edgelist.md)
or by `control$truncate.el.cuml`.

## See also

[`get_cumulative_degree()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_degree.md)
for a partner count per index;
[`get_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_edgelist.md)
/
[`get_cumulative_edgelists_df()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_edgelists_df.md)
for the underlying edgelist;
[`get_unique_ids()`](https://epimodel.github.io/EpiModel/reference/unique_id-tools.md)
/
[`get_posit_ids()`](https://epimodel.github.io/EpiModel/reference/unique_id-tools.md)
for ID conversion.

Other cumulative_edgelist:
[`as_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/as_cumulative_edgelist.md),
[`dedup_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/dedup_cumulative_edgelist.md),
[`get_cumulative_degree()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_degree.md),
[`get_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_edgelist.md),
[`get_cumulative_edgelists_df()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_edgelists_df.md),
[`reachable-nodes`](https://epimodel.github.io/EpiModel/reference/reachable-nodes.md),
[`update_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/update_cumulative_edgelist.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Contacts of nodes 1..10 across all networks within the last 30 steps,
# excluding partners who have since departed:
get_partners(dat,
             index_posit_ids = 1:10,
             truncate = 30,
             only.active.nodes = TRUE)
} # }
```
