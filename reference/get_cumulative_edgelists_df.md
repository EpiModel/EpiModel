# Get the Cumulative Edgelists of a Model

Get the Cumulative Edgelists of a Model

## Usage

``` r
get_cumulative_edgelists_df(dat, networks = NULL)
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

A `data.frame` with 5 columns:

- `index`: the unique ID (see `get_unique_ids`) of the indexes.

- `partner`: the unique ID (see `get_unique_ids`) of the
  partners/contacts.

- `start`: the time step in which the edge started.

- `stop`: the time step in which the edge stopped; if ongoing, then `NA`
  is returned.

- `network`: the numerical index for the network on which the
  partnership/contact is located.
