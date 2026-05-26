# Update a Cumulative Edgelist of the Specified Network

Records the current network state into the cumulative edgelist, stamping
newly observed edges with their `start` step and edges that have just
dissolved with their `stop` step. Returns `dat` unchanged when
cumulative-edgelist tracking is disabled.

## Usage

``` r
update_cumulative_edgelist(dat, network, truncate = 0)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md).

- network:

  Numerical index of the network for which the cumulative edgelist will
  be updated. (May be \> 1 for models with multiple overlapping
  networks.)

- truncate:

  After how many time steps a partnership that is no longer active
  should be removed from the output. See the Truncation section.

## Value

The updated `netsim_dat` main list object.

## Details

Calling this function is a no-op unless
`control.net(cumulative.edgelist = TRUE)`. With tracking enabled, the
built-in network-resimulation module
[`resim_nets()`](https://epimodel.github.io/EpiModel/reference/resim_nets.md)
calls this once per network at every time step, using
`truncate = control$truncate.el.cuml`. Custom modules that mutate the
network outside the TERGM machinery (e.g., forming or dissolving edges
directly) should call this manually after the change so the cumulative
record stays consistent.

## Truncation

To avoid storing a cumulative edgelist too long, the `truncate`
parameter defines a number of steps after which an edge that is no
longer active is dropped from the stored history. When `truncate = Inf`,
no edges are ever removed. When `truncate = 0` (the default), only
currently active edges are kept; this is useful for tracking the `start`
step of each active edge while keeping memory low.

## See also

[`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md)
(`cumulative.edgelist`, `truncate.el.cuml`),
[`get_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_edgelist.md)
to read the result.

Other cumulative_edgelist:
[`as_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/as_cumulative_edgelist.md),
[`dedup_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/dedup_cumulative_edgelist.md),
[`get_cumulative_degree()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_degree.md),
[`get_cumulative_edgelist()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_edgelist.md),
[`get_cumulative_edgelists_df()`](https://epimodel.github.io/EpiModel/reference/get_cumulative_edgelists_df.md),
[`get_partners()`](https://epimodel.github.io/EpiModel/reference/get_partners.md),
[`reachable-nodes`](https://epimodel.github.io/EpiModel/reference/reachable-nodes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Inside a custom module, after editing the network outside the TERGM:
for (n in seq_len(dat$num.nw)) {
  dat <- update_cumulative_edgelist(dat, n,
    truncate = get_control(dat, "truncate.el.cuml"))
}
} # }
```
