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
`(head, tail)` pair, sorted by pair and `start`.

## Details

Spells are inclusive on both ends, so two spells of a pair are merged
when one starts at or before the step the other stops; spells separated
by a gap are kept apart. A `stop` of `NA` marks a spell that is still
ongoing, and a merged spell that includes one stays ongoing. Columns
other than `start` and `stop` are treated as part of the pair, so spells
of the same pair on different layers of a `network` column are not
merged.

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
el <- data.frame(head = 1, tail = 2, start = c(1, 3, 10), stop = c(5, 7, NA))
dedup_cumulative_edgelist(el)
#>   head tail start stop
#> 1    1    2     1    7
#> 2    1    2    10   NA
```
