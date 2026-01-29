# Get a Cumulative Edgelist From a Specified Network

Get a Cumulative Edgelist From a Specified Network

## Usage

``` r
get_cumulative_edgelist(dat, network)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

- network:

  Numerical index of the network from which the cumulative edgelist
  should be extracted. (May be \> 1 for models with multiple overlapping
  networks.)

## Value

A cumulative edgelist in `data.frame` form with 4 columns:

- `head`: the unique ID (see `get_unique_ids`) of the head node on the
  edge.

- `tail`: the unique ID (see `get_unique_ids`) of the tail node on the
  edge.

- `start`: the time step in which the edge started.

- `stop`: the time step in which the edge stopped; if ongoing, then `NA`
  is returned.
