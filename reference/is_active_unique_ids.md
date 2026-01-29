# Are These Nodes Active (Unique IDs)

Are These Nodes Active (Unique IDs)

## Usage

``` r
is_active_unique_ids(dat, unique_ids)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

- unique_ids:

  A vector of node unique identifiers.

## Value

A logical vector with TRUE if the node is still active and FALSE
otherwise.
