# Are These Nodes Active (Positional IDs)

Are These Nodes Active (Positional IDs)

## Usage

``` r
is_active_posit_ids(dat, posit_ids)
```

## Arguments

- dat:

  Main `netsim_dat` object containing a `networkDynamic` object and
  other initialization information passed from
  [`netsim()`](http://epimodel.github.io/EpiModel/reference/netsim.md).

- posit_ids:

  A vector of node positional identifiers.

## Value

A logical vector with TRUE if the node is still active and FALSE
otherwise.
