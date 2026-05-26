# Deduplicate a Cumulative Edgelist

Merges overlapping spells for the same `(head, tail)` pair into a single
row. Useful when concatenating cumulative edgelists across multiple
simulations or sources, where the same partnership can appear in more
than one segment. The reachability functions assume non-overlapping
input.

## Usage

``` r
dedup_cumulative_edgelist(el)
```

## Arguments

- el:

  A cumulative edgelist with potentially overlapping edges.

## Value

A cumulative edgelist with no overlapping edges for the same
`(head, tail)` pair.

## See also

[`as_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/as_cumulative_edgelist.md)
for one source of duplicates;
[`get_forward_reachable()`](https://epimodel.github.io/EpiModel/reference/reachable-nodes.md)
/
[`get_backward_reachable()`](https://epimodel.github.io/EpiModel/reference/reachable-nodes.md),
which expect non-overlapping input.

Other cumulative_edgelist:
[`as_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/as_cumulative_edgelist.md),
[`get_cumulative_degree()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_degree.md),
[`get_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_edgelist.md),
[`get_cumulative_edgelists_df()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_edgelists_df.md),
[`get_partners()`](https://epimodel.github.io/EpiModel/reference/get_partners.md),
[`reachable-nodes`](https://epimodel.github.io/EpiModel/reference/reachable-nodes.md),
[`update_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/update_cumulative_edgelist.md)

## Examples

``` r
if (FALSE) { # \dontrun{
el_cuml <- dedup_cumulative_edgelist(el_raw)
} # }
```
