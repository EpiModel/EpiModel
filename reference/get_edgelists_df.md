# Get the Edgelist(s) from the Specified Network(s)

Get the Edgelist(s) from the Specified Network(s)

## Usage

``` r
get_edgelists_df(dat, networks = NULL)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

- networks:

  Numerical indexes of the networks to extract the partnerships from.
  (May be \> 1 for models with multiple overlapping networks.) If
  `NULL`, extract from all networks.

## Value

A `data.frame` with the following columns:

- `head`: Positional ID of the head node.

- `tail`: Positional ID of the tail node.

- `network`: The numerical index of the network on which the edge is
  located.
