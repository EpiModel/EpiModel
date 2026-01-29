# Update a Cumulative Edgelist of the Specified Network

Update a Cumulative Edgelist of the Specified Network

## Usage

``` r
update_cumulative_edgelist(dat, network, truncate = 0)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

- network:

  Numerical index of the network for which the cumulative edgelist will
  be updated. (May be \> 1 for models with multiple overlapping
  networks.)

- truncate:

  After how many time steps a partnership that is no longer active
  should be removed from the output.

## Value

The updated `netsim_dat` main list object.

## Truncation

To avoid storing a cumulative edgelist too long, the `truncate`
parameter defines a number of steps after which an edge that is no
longer active is truncated out of the cumulative edgelist. When
`truncate = Inf`, no edges are ever removed. When `truncate = 0`, only
the active edges are kept. You may want this behavior to keep track of
the active edges' start step.
