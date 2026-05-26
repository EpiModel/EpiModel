# Convert an Object to a Cumulative Edgelist

Converts an object that carries timed-edge information into a
`cumulative_edgelist`: a `data.frame` with `head`, `tail`, `start`,
`stop` columns matching the format used elsewhere in EpiModel. The
typical use is feeding a `networkDynamic` object (from a non-`tergmLite`
simulation, an external source, or a post-hoc construction) into the
reachability tools.

## Usage

``` r
as_cumulative_edgelist(x)
```

## Arguments

- x:

  An object to be converted to a cumulative edgelist.

## Value

A `cumulative_edgelist` object: a `data.frame` with at least the columns
`head`, `tail`, `start`, `stop`. Edges are active from time `start` to
time `stop` inclusive; `NA` in `stop` means the edge was never observed
to dissolve.

## Details

The `networkDynamic` method subtracts 1 from each spell's `terminus` so
that the resulting `start`/`stop` range is inclusive on both ends,
matching the convention used by
[`get_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_edgelist.md)
and the reachability functions.

## See also

[`get_forward_reachable()`](https://epimodel.github.io/EpiModel/reference/reachable-nodes.md)
/
[`get_backward_reachable()`](https://epimodel.github.io/EpiModel/reference/reachable-nodes.md)
for typical downstream uses;
[`dedup_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/dedup_cumulative_edgelist.md)
when the source contains overlapping spells for the same pair of nodes.

Other cumulative_edgelist:
[`dedup_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/dedup_cumulative_edgelist.md),
[`get_cumulative_degree()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_degree.md),
[`get_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_edgelist.md),
[`get_cumulative_edgelists_df()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_edgelists_df.md),
[`get_partners()`](https://epimodel.github.io/EpiModel/reference/get_partners.md),
[`reachable-nodes`](https://epimodel.github.io/EpiModel/reference/reachable-nodes.md),
[`update_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/update_cumulative_edgelist.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Convert a saved networkDynamic from a non-tergmLite netsim run
nd <- get_network(sim, sim = 1, network = 1)
el_cuml <- as_cumulative_edgelist(nd)
fwd <- get_forward_reachable(el_cuml, from_step = 1, to_step = 100)
} # }
```
